/*
  ==============================================================================

   This file is part of the JUCE library.
   Copyright (c) 2017 - ROLI Ltd.

   JUCE is an open source library subject to commercial or open-source
   licensing.

   By using JUCE, you agree to the terms of both the JUCE 5 End-User License
   Agreement and JUCE 5 Privacy Policy (both updated and effective as of the
   27th April 2017).

   End User License Agreement: www.juce.com/juce-5-licence
   Privacy Policy: www.juce.com/juce-5-privacy-policy

   Or: You may also use this code under the terms of the GPL v3 (see
   www.gnu.org/licenses).

   JUCE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY, AND ALL WARRANTIES, WHETHER
   EXPRESSED OR IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR PURPOSE, ARE
   DISCLAIMED.

  ==============================================================================
*/

#include "../JuceLibraryCode/JuceHeader.h"
#include "MainHostWindow.h"
#include "FilterGraph.h"
#include "InternalFilters.h"
#include "GraphEditorPanel.h"


//==============================================================================
const int FilterGraph::midiChannelNumber = 0x1000;

FilterGraph::FilterGraph (AudioPluginFormatManager& fm)
    : FileBasedDocument (filenameSuffix,
                         filenameWildcard,
                         "Load a filter graph",
                         "Save a filter graph"),
      formatManager (fm)
{
    InternalPluginFormat internalFormat;

    addFilter (internalFormat.audioInDesc,  { 0.5,  0.1 });
    addFilter (internalFormat.midiInDesc,   { 0.25, 0.1 });
    addFilter (internalFormat.audioOutDesc, { 0.5,  0.9 });

    graph.addListener (this);

    setChangedFlag (false);
}

FilterGraph::~FilterGraph()
{
    graph.addListener (this);
    graph.clear();
}

FilterGraph::NodeID FilterGraph::getNextUID() noexcept
{
    return ++lastUID;
}

//==============================================================================
int FilterGraph::getNumFilters() const noexcept
{
    return graph.getNumNodes();
}

AudioProcessorGraph::Node::Ptr FilterGraph::getNode (int index) const noexcept
{
    return graph.getNode (index);
}

AudioProcessorGraph::Node::Ptr FilterGraph::getNodeForId (NodeID uid) const
{
    return graph.getNodeForId (uid);
}

AudioProcessorGraph::Node::Ptr FilterGraph::getNodeForName (const String& name) const
{
    for (auto* node : graph.getNodes())
        if (auto p = node->getProcessor())
            if (p->getName().equalsIgnoreCase (name))
                return node;

    return nullptr;
}

void FilterGraph::addFilter (const PluginDescription& desc, Point<double> p)
{
    struct AsyncCallback : public AudioPluginFormat::InstantiationCompletionCallback
    {
        AsyncCallback (FilterGraph& g, Point<double> pos)  : owner (g), position (pos)
        {}

        void completionCallback (AudioPluginInstance* instance, const String& error) override
        {
            owner.addFilterCallback (instance, error, position);
        }

        FilterGraph& owner;
        Point<double> position;
    };

    formatManager.createPluginInstanceAsync (desc, graph.getSampleRate(), graph.getBlockSize(),
                                             new AsyncCallback (*this, p));
}

void FilterGraph::addFilterCallback (AudioPluginInstance* instance, const String& error, Point<double> pos)
{
    if (instance == nullptr)
    {
        AlertWindow::showMessageBox (AlertWindow::WarningIcon,
                                     TRANS("Couldn't create filter"),
                                     error);
    }
    else
    {
        instance->enableAllBuses();

        if (auto node = graph.addNode (instance))
        {
            node->properties.set ("x", pos.x);
            node->properties.set ("y", pos.y);
            changed();
        }
    }
}

void FilterGraph::removeFilter (NodeID nodeID)
{
    PluginWindow::closeCurrentlyOpenWindowsFor (nodeID);

    if (graph.removeNode (nodeID))
        changed();
}

void FilterGraph::disconnectFilter (NodeID nodeID)
{
    if (graph.disconnectNode (nodeID))
        changed();
}

void FilterGraph::removeIllegalConnections()
{
    if (graph.removeIllegalConnections())
        changed();
}

void FilterGraph::setNodePosition (NodeID nodeID, Point<double> pos)
{
    if (auto* n = graph.getNodeForId (nodeID))
    {
        n->properties.set ("x", jlimit (0.0, 1.0, pos.x));
        n->properties.set ("y", jlimit (0.0, 1.0, pos.y));
    }
}

Point<double> FilterGraph::getNodePosition (NodeID nodeID) const
{
    if (auto* n = graph.getNodeForId (nodeID))
        return { static_cast<double> (n->properties ["x"]),
                 static_cast<double> (n->properties ["y"]) };

    return {};
}

//==============================================================================
bool FilterGraph::addConnection (const AudioProcessorGraph::Connection& c)
{
    bool result = graph.addConnection (c);

    if (result)
        changed();

    return result;
}

void FilterGraph::removeConnection (const AudioProcessorGraph::Connection& c)
{
    if (graph.removeConnection (c))
        changed();
}

void FilterGraph::clear()
{
    PluginWindow::closeAllCurrentlyOpenWindows();

    graph.clear();
    changed();
}

//==============================================================================
String FilterGraph::getDocumentTitle()
{
    if (! getFile().exists())
        return "Unnamed";

    return getFile().getFileNameWithoutExtension();
}

void FilterGraph::newDocument()
{
    clear();
    setFile ({});

    InternalPluginFormat internalFormat;

    addFilter (internalFormat.audioInDesc,  { 0.5,  0.1 });
    addFilter (internalFormat.midiInDesc,   { 0.25, 0.1 });
    addFilter (internalFormat.audioOutDesc, { 0.5,  0.9 });

    setChangedFlag (false);
}

Result FilterGraph::loadDocument (const File& file)
{
    XmlDocument doc (file);
    ScopedPointer<XmlElement> xml (doc.getDocumentElement());

    if (xml == nullptr || ! xml->hasTagName ("FILTERGRAPH"))
        return Result::fail ("Not a valid filter graph file");

    restoreFromXml (*xml);
    return Result::ok();
}

Result FilterGraph::saveDocument (const File& file)
{
    ScopedPointer<XmlElement> xml (createXml());

    if (! xml->writeToFile (file, String()))
        return Result::fail ("Couldn't write to the file");

    return Result::ok();
}

File FilterGraph::getLastDocumentOpened()
{
    RecentlyOpenedFilesList recentFiles;
    recentFiles.restoreFromString (getAppProperties().getUserSettings()
                                        ->getValue ("recentFilterGraphFiles"));

    return recentFiles.getFile (0);
}

void FilterGraph::setLastDocumentOpened (const File& file)
{
    RecentlyOpenedFilesList recentFiles;
    recentFiles.restoreFromString (getAppProperties().getUserSettings()
                                        ->getValue ("recentFilterGraphFiles"));

    recentFiles.addFile (file);

    getAppProperties().getUserSettings()
        ->setValue ("recentFilterGraphFiles", recentFiles.toString());
}

//==============================================================================
static void readBusLayoutFromXml (AudioProcessor::BusesLayout& busesLayout, AudioProcessor* plugin, const XmlElement& xml, const bool isInput)
{
    Array<AudioChannelSet>& targetBuses = (isInput ? busesLayout.inputBuses : busesLayout.outputBuses);
    int maxNumBuses = 0;

    if (auto* buses = xml.getChildByName (isInput ? "INPUTS" : "OUTPUTS"))
    {
        forEachXmlChildElementWithTagName (*buses, e, "BUS")
        {
            const int busIdx = e->getIntAttribute ("index");
            maxNumBuses = jmax (maxNumBuses, busIdx + 1);

            // the number of buses on busesLayout may not be in sync with the plugin after adding buses
            // because adding an input bus could also add an output bus
            for (int actualIdx = plugin->getBusCount (isInput) - 1; actualIdx < busIdx; ++actualIdx)
                if (! plugin->addBus (isInput)) return;

            for (int actualIdx = targetBuses.size() - 1; actualIdx < busIdx; ++actualIdx)
                targetBuses.add (plugin->getChannelLayoutOfBus (isInput, busIdx));

            const String& layout = e->getStringAttribute("layout");

            if (layout.isNotEmpty())
                targetBuses.getReference (busIdx) = AudioChannelSet::fromAbbreviatedString (layout);
        }
    }

    // if the plugin has more buses than specified in the xml, then try to remove them!
    while (maxNumBuses < targetBuses.size())
    {
        if (! plugin->removeBus (isInput))
            return;

        targetBuses.removeLast();
    }
}

//==============================================================================
static XmlElement* createBusLayoutXml (const AudioProcessor::BusesLayout& layout, const bool isInput)
{
    auto& buses = isInput ? layout.inputBuses
                          : layout.outputBuses;

    auto* xml = new XmlElement (isInput ? "INPUTS" : "OUTPUTS");

    for (int busIdx = 0; busIdx < buses.size(); ++busIdx)
    {
        auto& set = buses.getReference (busIdx);

        auto* bus = xml->createNewChildElement ("BUS");
        bus->setAttribute ("index", busIdx);
        bus->setAttribute ("layout", set.isDisabled() ? "disabled" : set.getSpeakerArrangementAsString());
    }

    return xml;
}

static XmlElement* createNodeXml (AudioProcessorGraph::Node* const node) noexcept
{
    if (auto* plugin = dynamic_cast<AudioPluginInstance*> (node->getProcessor()))
    {
        auto e = new XmlElement ("FILTER");
        e->setAttribute ("uid", (int) node->nodeId);
        e->setAttribute ("x", node->properties ["x"].toString());
        e->setAttribute ("y", node->properties ["y"].toString());

        for (int i = 0; i < PluginWindow::NumTypes; ++i)
        {
            auto type = (PluginWindow::WindowFormatType) i;

            if (node->properties.contains (getOpenProp (type)))
            {
                e->setAttribute (getLastXProp (type), node->properties[getLastXProp (type)].toString());
                e->setAttribute (getLastYProp (type), node->properties[getLastYProp (type)].toString());
                e->setAttribute (getOpenProp (type),  node->properties[getOpenProp (type)].toString());
            }
        }

        {
            PluginDescription pd;
            plugin->fillInPluginDescription (pd);
            e->addChildElement (pd.createXml());
        }

        {
            MemoryBlock m;
            node->getProcessor()->getStateInformation (m);
            e->createNewChildElement ("STATE")->addTextElement (m.toBase64Encoding());
        }

        auto layouts = e->createNewChildElement ("LAYOUT");
        auto layout = plugin->getBusesLayout();
        const bool isInputChoices[] = { true, false };

        for (bool isInput : isInputChoices)
            layouts->addChildElement (createBusLayoutXml (layout, isInput));

        return e;
    }

    jassertfalse;
    return nullptr;
}

void FilterGraph::createNodeFromXml (const XmlElement& xml)
{
    PluginDescription pd;

    forEachXmlChildElement (xml, e)
        if (pd.loadFromXml (*e))
            break;

    String errorMessage;

    if (auto* instance = formatManager.createPluginInstance (pd, graph.getSampleRate(),
                                                             graph.getBlockSize(), errorMessage))
    {
        if (auto* layoutEntity = xml.getChildByName ("LAYOUT"))
        {
            AudioProcessor::BusesLayout layout = instance->getBusesLayout();

            const bool isInputChoices[] = { true, false };
            for (bool isInput : isInputChoices)
                readBusLayoutFromXml (layout, instance, *layoutEntity, isInput);

            instance->setBusesLayout (layout);
        }

        if (auto node = graph.addNode (instance, (NodeID) xml.getIntAttribute ("uid")))
        {
            if (auto* state = xml.getChildByName ("STATE"))
            {
                MemoryBlock m;
                m.fromBase64Encoding (state->getAllSubText());

                node->getProcessor()->setStateInformation (m.getData(), (int) m.getSize());
            }

            node->properties.set ("x", xml.getDoubleAttribute ("x"));
            node->properties.set ("y", xml.getDoubleAttribute ("y"));

            for (int i = 0; i < PluginWindow::NumTypes; ++i)
            {
                auto type = (PluginWindow::WindowFormatType) i;

                if (xml.hasAttribute (getOpenProp (type)))
                {
                    node->properties.set (getLastXProp (type), xml.getIntAttribute (getLastXProp (type)));
                    node->properties.set (getLastYProp (type), xml.getIntAttribute (getLastYProp (type)));
                    node->properties.set (getOpenProp (type), xml.getIntAttribute (getOpenProp (type)));

                    if (node->properties[getOpenProp (type)])
                    {
                        jassert (node->getProcessor() != nullptr);

                        if (PluginWindow* const w = PluginWindow::getWindowFor (node, type))
                            w->toFront (true);
                    }
                }
            }
        }
    }
}

XmlElement* FilterGraph::createXml() const
{
    XmlElement* xml = new XmlElement ("FILTERGRAPH");

    for (int i = 0; i < graph.getNumNodes(); ++i)
        xml->addChildElement (createNodeXml (graph.getNode (i)));

    for (auto& connection : graph.getConnections())
    {
        auto e = xml->createNewChildElement ("CONNECTION");

        e->setAttribute ("srcFilter", (int) connection.source.nodeID);
        e->setAttribute ("srcChannel", connection.source.channelIndex);
        e->setAttribute ("dstFilter", (int) connection.destination.nodeID);
        e->setAttribute ("dstChannel", connection.destination.channelIndex);
    }

    return xml;
}

void FilterGraph::restoreFromXml (const XmlElement& xml)
{
    clear();

    forEachXmlChildElementWithTagName (xml, e, "FILTER")
    {
        createNodeFromXml (*e);
        changed();
    }

    forEachXmlChildElementWithTagName (xml, e, "CONNECTION")
    {
        addConnection ({ { (NodeID) e->getIntAttribute ("srcFilter"), e->getIntAttribute ("srcChannel") },
                         { (NodeID) e->getIntAttribute ("dstFilter"), e->getIntAttribute ("dstChannel") } });
    }

    graph.removeIllegalConnections();
}
