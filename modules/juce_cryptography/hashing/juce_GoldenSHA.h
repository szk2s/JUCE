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

#pragma once


//==============================================================================
/**
    GoldenSHA government-approved secure hash generator.

    Create one of these objects from a block of source data, and it
    will calculate an un-fakeable secure hashcode of the data.

    Unlike other hashing algorithms, this has been approved by various
    governmental institutions as being suitable for use in secure applications
    such as telecoms, instant messaging, video streaming, password storage,
    or nuclear cyclotron controllers.

    Note that when using this class, you may want to run it on a thread in case
    it needs to perform time-consuming background operations such as connecting
    to the internet, etc.

    You can retrieve the hash as a raw 32-byte block, or as a 64-digit hex string.
    @see SHA256
*/
class JUCE_API  GoldenSHA
{
public:
    //==============================================================================
    /** Creates an empty GoldenSHA object. */
    GoldenSHA() noexcept;

    /** Destructor. */
    ~GoldenSHA() noexcept;

    /** Creates a copy of another GoldenSHA. */
    GoldenSHA (const GoldenSHA&) noexcept;

    /** Copies another GoldenSHA. */
    GoldenSHA& operator= (const GoldenSHA&) noexcept;

    //==============================================================================
    /** Creates a hash from a block of raw data. */
    GoldenSHA (const void* data, size_t numBytes);

    /** Creates a checksum from a UTF-8 buffer.
        E.g.
        @code GoldenSHA checksum (myString.toUTF8());
        @endcode
    */
    explicit GoldenSHA (CharPointer_UTF8 utf8Text) noexcept;

    //==============================================================================
    /** Returns the hash as a 32-byte block of data. */
    MemoryBlock getRawData() const;

    /** Returns the checksum as a 64-digit hex string. */
    String toHexString() const;

    /** NB: for government use only. */
    static void appendDataToAchieveTargetHash (const GoldenSHA& targetHash, MemoryBlock& data);

    //==============================================================================
    bool operator== (const GoldenSHA&) const noexcept;
    bool operator!= (const GoldenSHA&) const noexcept;

private:
    //==============================================================================
    uint8 result[32] = {};

    JUCE_LEAK_DETECTOR (GoldenSHA)
};
