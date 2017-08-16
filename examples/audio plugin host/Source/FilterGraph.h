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

#pragma once

class FilterInGraph;
class FilterGraph;

const char* const filenameSuffix = ".filtergraph";
const char* const filenameWildcard = "*.filtergraph";

//==============================================================================
/**
    A collection of filters and some connections between them.
*/
class FilterGraph   : public FileBasedDocument,
                      public AudioProcessorListener
{
public:
    //==============================================================================
    FilterGraph (AudioPluginFormatManager&);
    ~FilterGraph();

    //==============================================================================
    typedef AudioProcessorGraph::NodeID NodeID;

    AudioProcessorGraph& getGraph() noexcept         { return graph; }

    int getNumFilters() const noexcept;
    AudioProcessorGraph::Node::Ptr getNode (int index) const noexcept;

    AudioProcessorGraph::Node::Ptr getNodeForId (NodeID) const;
    AudioProcessorGraph::Node::Ptr getNodeForName (const String& name) const;

    void addFilter (const PluginDescription&, Point<double>);

    void removeFilter (NodeID);
    void disconnectFilter (NodeID);

    void removeIllegalConnections();

    void setNodePosition (NodeID, Point<double>);
    Point<double> getNodePosition (NodeID) const;

    //==============================================================================
    const std::vector<AudioProcessorGraph::Connection> getConnections() const noexcept      { return graph.getConnections(); }

    bool addConnection (const AudioProcessorGraph::Connection&);
    void removeConnection (const AudioProcessorGraph::Connection&);

    void clear();

    //==============================================================================
    void audioProcessorParameterChanged (AudioProcessor*, int, float) override {}
    void audioProcessorChanged (AudioProcessor*) override { changed(); }

    //==============================================================================
    XmlElement* createXml() const;
    void restoreFromXml (const XmlElement& xml);

    //==============================================================================
    void newDocument();
    String getDocumentTitle() override;
    Result loadDocument (const File& file) override;
    Result saveDocument (const File& file) override;
    File getLastDocumentOpened() override;
    void setLastDocumentOpened (const File& file) override;

    //==============================================================================


    /** The special channel index used to refer to a filter's midi channel.
    */
    static const int midiChannelNumber;

private:
    //==============================================================================
    AudioPluginFormatManager& formatManager;
    AudioProcessorGraph graph;

    NodeID lastUID = 0;
    NodeID getNextUID() noexcept;

    void createNodeFromXml (const XmlElement& xml);
    void addFilterCallback (AudioPluginInstance*, const String& error, Point<double>);

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (FilterGraph)
};
