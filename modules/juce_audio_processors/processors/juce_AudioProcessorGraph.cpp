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


template <typename FloatType>
struct GraphRenderSequence
{
    GraphRenderSequence() {}

    struct Context
    {
        FloatType** audioBuffers;
        MidiBuffer* midiBuffers;
        int numSamples;
    };

    void perform (AudioBuffer<FloatType>& buffer, MidiBuffer& midiMessages)
    {
        auto numSamples = buffer.getNumSamples();
        auto maxSamples = renderingBuffer.getNumSamples();

        if (numSamples > maxSamples)
        {
            // being asked to render more samples than our buffers have, so slice things up...
            tempMIDI.clear();
            tempMIDI.addEvents (midiMessages, maxSamples, numSamples, -maxSamples);

            {
                AudioBuffer<FloatType> startAudio (buffer.getArrayOfWritePointers(), buffer.getNumChannels(), maxSamples);
                midiMessages.clear (maxSamples, numSamples);
                perform (startAudio, midiMessages);
            }

            AudioBuffer<FloatType> endAudio (buffer.getArrayOfWritePointers(), buffer.getNumChannels(), maxSamples, numSamples - maxSamples);
            perform (endAudio, tempMIDI);
            return;
        }

        currentAudioInputBuffer = &buffer;
        currentAudioOutputBuffer.setSize (jmax (1, buffer.getNumChannels()), numSamples);
        currentAudioOutputBuffer.clear();
        currentMidiInputBuffer = &midiMessages;
        currentMidiOutputBuffer.clear();

        {
            const Context context { renderingBuffer.getArrayOfWritePointers(), midiBuffers.begin(), numSamples };

            for (auto* op : renderOps)
                op->perform (context);
        }

        for (int i = 0; i < buffer.getNumChannels(); ++i)
            buffer.copyFrom (i, 0, currentAudioOutputBuffer, i, 0, numSamples);

        midiMessages.clear();
        midiMessages.addEvents (currentMidiOutputBuffer, 0, buffer.getNumSamples(), 0);
        currentAudioInputBuffer = nullptr;
    }

    void addClearChannelOp (int index)
    {
        createOp ([=] (const Context& c)    { FloatVectorOperations::clear (c.audioBuffers[index], c.numSamples); });
    }

    void addCopyChannelOp (int srcIndex, int dstIndex)
    {
        createOp ([=] (const Context& c)    { FloatVectorOperations::copy (c.audioBuffers[dstIndex],
                                                                           c.audioBuffers[srcIndex],
                                                                           c.numSamples); });
    }

    void addAddChannelOp (int srcIndex, int dstIndex)
    {
        createOp ([=] (const Context& c)    { FloatVectorOperations::add (c.audioBuffers[dstIndex],
                                                                          c.audioBuffers[srcIndex],
                                                                          c.numSamples); });
    }

    void addClearMidiBufferOp (int index)
    {
        createOp ([=] (const Context& c)    { c.midiBuffers[index].clear(); });
    }

    void addCopyMidiBufferOp (int srcIndex, int dstIndex)
    {
        createOp ([=] (const Context& c)    { c.midiBuffers[dstIndex] = c.midiBuffers[srcIndex]; });
    }

    void addAddMidiBufferOp (int srcIndex, int dstIndex)
    {
        createOp ([=] (const Context& c)    { c.midiBuffers[dstIndex].addEvents (c.midiBuffers[srcIndex],
                                                                                 0, c.numSamples, 0); });
    }

    void addDelayChannelOp (int chan, int delaySize)
    {
        renderOps.add (new DelayChannelOp (chan, delaySize));
    }

    void addProcessOp (const AudioProcessorGraph::Node::Ptr& node,
                       const Array<int>& audioChannelsUsed, int totalNumChans, int midiBuffer)
    {
        renderOps.add (new ProcessOp (node, audioChannelsUsed, totalNumChans, midiBuffer));
    }

    void prepareBuffers (int blockSize)
    {
        renderingBuffer.setSize (numBuffersNeeded + 1, blockSize);
        renderingBuffer.clear();
        currentAudioOutputBuffer.setSize (numBuffersNeeded + 1, blockSize);
        currentAudioOutputBuffer.clear();

        currentAudioInputBuffer = nullptr;
        currentMidiInputBuffer = nullptr;
        currentMidiOutputBuffer.clear();

        midiBuffers.clearQuick();
        midiBuffers.resize (numMidiBuffersNeeded);

        const int defaultMIDIBufferSize = 512;

        tempMIDI.ensureSize (defaultMIDIBufferSize);

        for (auto&& m : midiBuffers)
            m.ensureSize (defaultMIDIBufferSize);
    }

    void releaseBuffers()
    {
        renderingBuffer.setSize (1, 1);
        currentAudioOutputBuffer.setSize (1, 1);
        currentAudioInputBuffer = nullptr;
        currentMidiInputBuffer = nullptr;
        currentMidiOutputBuffer.clear();
        midiBuffers.clear();
    }

    int numBuffersNeeded = 0, numMidiBuffersNeeded = 0;

    AudioBuffer<FloatType> renderingBuffer, currentAudioOutputBuffer;
    AudioBuffer<FloatType>* currentAudioInputBuffer = nullptr;

    MidiBuffer* currentMidiInputBuffer = nullptr;
    MidiBuffer currentMidiOutputBuffer;

    Array<MidiBuffer> midiBuffers;
    MidiBuffer tempMIDI;

private:
    //==============================================================================
    struct RenderingOp
    {
        RenderingOp() noexcept {}
        virtual ~RenderingOp() {}
        virtual void perform (const Context&) = 0;

        JUCE_LEAK_DETECTOR (RenderingOp)
    };

    OwnedArray<RenderingOp> renderOps;

    //==============================================================================
    template <typename LambdaType>
    void createOp (LambdaType&& fn)
    {
        struct LambdaOp  : public RenderingOp
        {
            LambdaOp (LambdaType&& f) : function (static_cast<LambdaType&&> (f)) {}
            void perform (const Context& c) override    { function (c); }

            LambdaType function;
        };

        renderOps.add (new LambdaOp (static_cast<LambdaType&&> (fn)));
    }

    //==============================================================================
    struct DelayChannelOp  : public RenderingOp
    {
        DelayChannelOp (int chan, int delaySize)
            : channel (chan),
              bufferSize (delaySize + 1),
              writeIndex (delaySize)
        {
            buffer.calloc ((size_t) bufferSize);
        }

        void perform (const Context& c)
        {
            auto* data = c.audioBuffers[channel];

            for (int i = c.numSamples; --i >= 0;)
            {
                buffer[writeIndex] = *data;
                *data++ = buffer[readIndex];

                if (++readIndex  >= bufferSize) readIndex = 0;
                if (++writeIndex >= bufferSize) writeIndex = 0;
            }
        }

        HeapBlock<FloatType> buffer;
        const int channel, bufferSize;
        int readIndex = 0, writeIndex;

        JUCE_DECLARE_NON_COPYABLE (DelayChannelOp)
    };

    //==============================================================================
    struct ProcessOp   : public RenderingOp
    {
        ProcessOp (const AudioProcessorGraph::Node::Ptr& n,
                   const Array<int>& audioChannelsUsed,
                   int totalNumChans, int midiBuffer)
            : node (n),
              processor (*n->getProcessor()),
              audioChannelsToUse (audioChannelsUsed),
              totalChans (jmax (1, totalNumChans)),
              midiBufferToUse (midiBuffer)
        {
            audioChannels.calloc ((size_t) totalChans);

            while (audioChannelsToUse.size() < totalChans)
                audioChannelsToUse.add (0);
        }

        void perform (const Context& c)
        {
            for (int i = 0; i < totalChans; ++i)
                audioChannels[i] = c.audioBuffers[audioChannelsToUse.getUnchecked (i)];

            AudioBuffer<FloatType> buffer (audioChannels, totalChans, c.numSamples);

            if (processor.isSuspended())
                buffer.clear();
            else
                callProcess (buffer, c.midiBuffers[midiBufferToUse]);
        }

        void callProcess (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
        {
            if (processor.isUsingDoublePrecision())
            {
                tempBufferDouble.makeCopyOf (buffer, true);
                processor.processBlock (tempBufferDouble, midiMessages);
                buffer.makeCopyOf (tempBufferDouble, true);
            }
            else
            {
                processor.processBlock (buffer, midiMessages);
            }
        }

        void callProcess (AudioBuffer<double>& buffer, MidiBuffer& midiMessages)
        {
            if (processor.isUsingDoublePrecision())
            {
                processor.processBlock (buffer, midiMessages);
            }
            else
            {
                tempBufferFloat.makeCopyOf (buffer, true);
                processor.processBlock (tempBufferFloat, midiMessages);
                buffer.makeCopyOf (tempBufferFloat, true);
            }
        }

        const AudioProcessorGraph::Node::Ptr node;
        AudioProcessor& processor;

        Array<int> audioChannelsToUse;
        HeapBlock<FloatType*> audioChannels;
        AudioBuffer<float> tempBufferFloat, tempBufferDouble;
        const int totalChans, midiBufferToUse;

        JUCE_DECLARE_NON_COPYABLE (ProcessOp)
    };
};

//==============================================================================
//==============================================================================
template <typename RenderSequence>
struct RenderSequenceBuilder
{
    RenderSequenceBuilder (AudioProcessorGraph& g, RenderSequence& s)
        : graph (g), sequence (s)
    {
        createOrderedNodeList();

        audioBuffers.add (AssignedBuffer::createReadOnlyEmpty()); // first buffer is read-only zeros
        midiBuffers .add (AssignedBuffer::createReadOnlyEmpty());

        for (int i = 0; i < orderedNodes.size(); ++i)
        {
            createRenderingOpsForNode (*orderedNodes.getUnchecked(i), i);
            markAnyUnusedBuffersAsFree (audioBuffers, i);
            markAnyUnusedBuffersAsFree (midiBuffers, i);
        }

        graph.setLatencySamples (totalLatency);

        s.numBuffersNeeded = audioBuffers.size();
        s.numMidiBuffersNeeded = midiBuffers.size();
    }

    //==============================================================================
    typedef AudioProcessorGraph::NodeID NodeID;

    AudioProcessorGraph& graph;
    RenderSequence& sequence;

    Array<AudioProcessorGraph::Node*> orderedNodes;

    struct AssignedBuffer
    {
        AudioProcessorGraph::NodeAndChannel channel;

        static AssignedBuffer createReadOnlyEmpty() noexcept    { return { { (NodeID) zeroNodeID, 0 } }; }
        static AssignedBuffer createFree() noexcept             { return { { (NodeID) freeNodeID, 0 } }; }

        bool isReadOnlyEmpty() const noexcept                   { return channel.nodeID == (NodeID) zeroNodeID; }
        bool isFree() const noexcept                            { return channel.nodeID == (NodeID) freeNodeID; }
        bool isAssigned() const noexcept                        { return ! (isReadOnlyEmpty() || isFree()); }

        void setFree() noexcept                                 { channel = { (NodeID) freeNodeID, 0 }; }
        void setAssignedToNonExistentNode() noexcept            { channel = { (NodeID) anonNodeID, 0 }; }

    private:
        enum
        {
            anonNodeID = 0x7ffffffd,
            zeroNodeID = 0x7ffffffe,
            freeNodeID = 0x7fffffff
        };
    };

    Array<AssignedBuffer> audioBuffers, midiBuffers;

    static int getReadOnlyEmptyBufferIndex() { return 0; }

    struct Delay
    {
        NodeID nodeID;
        int delay;
    };

    std::unordered_map<NodeID, int> delays;
    int totalLatency = 0;

    int getNodeDelay (NodeID nodeID) const
    {
        auto n = delays.find (nodeID);
        return n != delays.end() ? n->second : 0;
    }

    int getInputLatencyForNode (NodeID nodeID) const
    {
        int maxLatency = 0;

        for (auto&& c : graph.getConnections())
            if (c.destination.nodeID == nodeID)
                maxLatency = jmax (maxLatency, getNodeDelay (c.source.nodeID));

        return maxLatency;
    }

    //==============================================================================
    // Fast lookup table for checking which nodes are inputs to others.
    struct NodeInputLookupTable
    {
        template <typename ConnectionList>
        explicit NodeInputLookupTable (const ConnectionList& connections)
        {
            for (auto&& c : connections)
                entries[c.destination.nodeID].add (c.source.nodeID);
        }

        bool isAnInputTo (NodeID sourceID, NodeID destID) const noexcept
        {
            return isAnInputToRecursive (sourceID, destID, entries.size());
        }

        bool isAnInputToRecursive (NodeID sourceID, NodeID destID, size_t recursionDepth) const noexcept
        {
            auto&& entry = entries.find (destID);

            if (entry != entries.end())
            {
                auto& srcNodes = entry->second;

                if (srcNodes.contains (sourceID))
                    return true;

                if (recursionDepth != 0)
                    for (auto&& src : srcNodes)
                        if (isAnInputToRecursive (sourceID, src, recursionDepth - 1))
                            return true;
            }

            return false;
        }

        std::unordered_map<NodeID, SortedSet<NodeID>> entries;
    };

    void createOrderedNodeList()
    {
        const NodeInputLookupTable table (graph.getConnections());

        for (auto* node : graph.getNodes())
        {
            int j = 0;

            for (; j < orderedNodes.size(); ++j)
                if (table.isAnInputTo (node->nodeId, orderedNodes.getUnchecked(j)->nodeId))
                  break;

            orderedNodes.insert (j, node);
        }
    }

    int findBufferForInputAudioChannel (AudioProcessorGraph::Node& node, const int inputChan,
                                        const int ourRenderingIndex, const int maxLatency)
    {
        auto& processor = *node.getProcessor();
        auto numOuts = processor.getTotalNumOutputChannels();

        auto sources = getSourcesForChannel (node, inputChan);

        // Handle an unconnected input channel...
        if (sources.isEmpty())
        {
            if (inputChan >= numOuts)
                return getReadOnlyEmptyBufferIndex();

            auto index = getFreeBuffer (audioBuffers);
            sequence.addClearChannelOp (index);
            return index;
        }

        // Handle an input from a single source..
        if (sources.size() == 1)
        {
            // channel with a straightforward single input..
            auto src = sources.getUnchecked(0);

            int bufIndex = getBufferContaining (src);

            if (bufIndex < 0)
            {
                // if not found, this is probably a feedback loop
                bufIndex = getReadOnlyEmptyBufferIndex();
                jassert (bufIndex >= 0);
            }

            if (inputChan < numOuts
                 && isBufferNeededLater (ourRenderingIndex, inputChan, src))
            {
                // can't mess up this channel because it's needed later by another node,
                // so we need to use a copy of it..
                auto newFreeBuffer = getFreeBuffer (audioBuffers);
                sequence.addCopyChannelOp (bufIndex, newFreeBuffer);
                bufIndex = newFreeBuffer;
            }

            auto nodeDelay = getNodeDelay (src.nodeID);

            if (nodeDelay < maxLatency)
                sequence.addDelayChannelOp (bufIndex, maxLatency - nodeDelay);

            return bufIndex;
        }

        // Handle a mix of several outputs coming into this input..
        int reusableInputIndex = -1;
        int bufIndex = -1;

        for (int i = 0; i < sources.size(); ++i)
        {
            auto src = sources.getReference(i);
            auto sourceBufIndex = getBufferContaining (src);

            if (sourceBufIndex >= 0 && ! isBufferNeededLater (ourRenderingIndex, inputChan, src))
            {
                // we've found one of our input chans that can be re-used..
                reusableInputIndex = i;
                bufIndex = sourceBufIndex;

                auto nodeDelay = getNodeDelay (src.nodeID);

                if (nodeDelay < maxLatency)
                    sequence.addDelayChannelOp (bufIndex, maxLatency - nodeDelay);

                break;
            }
        }

        if (reusableInputIndex < 0)
        {
            // can't re-use any of our input chans, so get a new one and copy everything into it..
            bufIndex = getFreeBuffer (audioBuffers);
            jassert (bufIndex != 0);

            audioBuffers.getReference (bufIndex).setAssignedToNonExistentNode();

            auto srcIndex = getBufferContaining (sources.getFirst());

            if (srcIndex < 0)
                sequence.addClearChannelOp (bufIndex);  // if not found, this is probably a feedback loop
            else
                sequence.addCopyChannelOp (srcIndex, bufIndex);

            reusableInputIndex = 0;
            auto nodeDelay = getNodeDelay (sources.getFirst().nodeID);

            if (nodeDelay < maxLatency)
                sequence.addDelayChannelOp (bufIndex, maxLatency - nodeDelay);
        }

        for (int i = 0; i < sources.size(); ++i)
        {
            if (i != reusableInputIndex)
            {
                auto src = sources.getReference(i);
                int srcIndex = getBufferContaining (src);

                if (srcIndex >= 0)
                {
                    auto nodeDelay = getNodeDelay (src.nodeID);

                    if (nodeDelay < maxLatency)
                    {
                        if (! isBufferNeededLater (ourRenderingIndex, inputChan, src))
                        {
                            sequence.addDelayChannelOp (srcIndex, maxLatency - nodeDelay);
                        }
                        else // buffer is reused elsewhere, can't be delayed
                        {
                            auto bufferToDelay = getFreeBuffer (audioBuffers);
                            sequence.addCopyChannelOp (srcIndex, bufferToDelay);
                            sequence.addDelayChannelOp (bufferToDelay, maxLatency - nodeDelay);
                            srcIndex = bufferToDelay;
                        }
                    }

                    sequence.addAddChannelOp (srcIndex, bufIndex);
                }
            }
        }

        return bufIndex;
    }

    int findBufferForInputMidiChannel (AudioProcessorGraph::Node& node, const int ourRenderingIndex, const int maxLatency)
    {
        auto& processor = *node.getProcessor();
        auto sources = getSourcesForChannel (node, AudioProcessorGraph::midiChannelIndex);

        // No midi inputs..
        if (sources.isEmpty())
        {
            auto midiBufferToUse = getFreeBuffer (midiBuffers); // need to pick a buffer even if the processor doesn't use midi

            if (processor.acceptsMidi() || processor.producesMidi())
                sequence.addClearMidiBufferOp (midiBufferToUse);

            return midiBufferToUse;
        }

        // One midi input..
        if (sources.size() == 1)
        {
            auto src = sources.getReference (0);
            auto midiBufferToUse = getBufferContaining (src);

            if (midiBufferToUse >= 0)
            {
                if (isBufferNeededLater (ourRenderingIndex, AudioProcessorGraph::midiChannelIndex, src))
                {
                    // can't mess up this channel because it's needed later by another node, so we
                    // need to use a copy of it..
                    auto newFreeBuffer = getFreeBuffer (midiBuffers);
                    sequence.addCopyMidiBufferOp (midiBufferToUse, newFreeBuffer);
                    midiBufferToUse = newFreeBuffer;
                }
            }
            else
            {
                // probably a feedback loop, so just use an empty one..
                midiBufferToUse = getFreeBuffer (midiBuffers); // need to pick a buffer even if the processor doesn't use midi
            }

            return midiBufferToUse;
        }

        // Multiple midi inputs..
        int midiBufferToUse = -1;
        int reusableInputIndex = -1;

        for (int i = 0; i < sources.size(); ++i)
        {
            auto src = sources.getReference (i);
            auto sourceBufIndex = getBufferContaining (src);

            if (sourceBufIndex >= 0
                 && ! isBufferNeededLater (ourRenderingIndex, AudioProcessorGraph::midiChannelIndex, src))
            {
                // we've found one of our input buffers that can be re-used..
                reusableInputIndex = i;
                midiBufferToUse = sourceBufIndex;
                break;
            }
        }

        if (reusableInputIndex < 0)
        {
            // can't re-use any of our input buffers, so get a new one and copy everything into it..
            midiBufferToUse = getFreeBuffer (midiBuffers);
            jassert (midiBufferToUse >= 0);

            auto srcIndex = getBufferContaining (sources.getUnchecked(0));

            if (srcIndex >= 0)
                sequence.addCopyMidiBufferOp (srcIndex, midiBufferToUse);
            else
                sequence.addClearMidiBufferOp (midiBufferToUse);

            reusableInputIndex = 0;
        }

        for (int i = 0; i < sources.size(); ++i)
        {
            if (i != reusableInputIndex)
            {
                auto srcIndex = getBufferContaining (sources.getUnchecked(i));

                if (srcIndex >= 0)
                    sequence.addAddMidiBufferOp (srcIndex, midiBufferToUse);
            }
        }

        return midiBufferToUse;
    }

    void createRenderingOpsForNode (AudioProcessorGraph::Node& node, const int ourRenderingIndex)
    {
        auto& processor = *node.getProcessor();
        auto numIns  = processor.getTotalNumInputChannels();
        auto numOuts = processor.getTotalNumOutputChannels();
        auto totalChans = jmax (numIns, numOuts);

        Array<int> audioChannelsToUse;
        auto maxLatency = getInputLatencyForNode (node.nodeId);

        for (int inputChan = 0; inputChan < numIns; ++inputChan)
        {
            // get a list of all the inputs to this node
            auto index = findBufferForInputAudioChannel (node, inputChan, ourRenderingIndex, maxLatency);
            jassert (index >= 0);

            audioChannelsToUse.add (index);

            if (inputChan < numOuts)
                audioBuffers.getReference (index).channel = { node.nodeId, inputChan };
        }

        for (int outputChan = numIns; outputChan < numOuts; ++outputChan)
        {
            auto index = getFreeBuffer (audioBuffers);
            jassert (index != 0);
            audioChannelsToUse.add (index);

            audioBuffers.getReference (index).channel = { node.nodeId, outputChan };
        }

        auto midiBufferToUse = findBufferForInputMidiChannel (node, ourRenderingIndex, maxLatency);

        if (processor.producesMidi())
            midiBuffers.getReference (midiBufferToUse).channel = { node.nodeId, AudioProcessorGraph::midiChannelIndex };

        delays[node.nodeId] = maxLatency + processor.getLatencySamples();

        if (numOuts == 0)
            totalLatency = maxLatency;

        sequence.addProcessOp (&node, audioChannelsToUse, totalChans, midiBufferToUse);
    }

    //==============================================================================
    Array<AudioProcessorGraph::NodeAndChannel> getSourcesForChannel (AudioProcessorGraph::Node& node, int inputChannelIndex)
    {
        Array<AudioProcessorGraph::NodeAndChannel> results;
        AudioProcessorGraph::NodeAndChannel nc { node.nodeId, inputChannelIndex };

        for (auto&& c : graph.getConnections())
            if (c.destination == nc)
                results.add (c.source);

        return results;
    }

    static int getFreeBuffer (Array<AssignedBuffer>& buffers)
    {
        for (int i = 1; i < buffers.size(); ++i)
            if (buffers.getReference(i).isFree())
                return i;

        buffers.add (AssignedBuffer::createFree());
        return buffers.size() - 1;
    }

    int getBufferContaining (AudioProcessorGraph::NodeAndChannel output) const noexcept
    {
        int i = 0;

        for (auto& b : output.isMIDI() ? midiBuffers : audioBuffers)
        {
            if (b.channel == output)
                return i;

            ++i;
        }

        return -1;
    }

    void markAnyUnusedBuffersAsFree (Array<AssignedBuffer>& buffers, const int stepIndex)
    {
        for (auto& b : buffers)
            if (b.isAssigned() && ! isBufferNeededLater (stepIndex, -1, b.channel))
                b.setFree();
    }

    bool isBufferNeededLater (int stepIndexToSearchFrom,
                              int inputChannelOfIndexToIgnore,
                              AudioProcessorGraph::NodeAndChannel output) const
    {
        while (stepIndexToSearchFrom < orderedNodes.size())
        {
            auto* node = orderedNodes.getUnchecked (stepIndexToSearchFrom);

            if (output.isMIDI())
            {
                if (inputChannelOfIndexToIgnore != AudioProcessorGraph::midiChannelIndex
                    && graph.isConnected ({ { output.nodeID, AudioProcessorGraph::midiChannelIndex },
                                            { node->nodeId,  AudioProcessorGraph::midiChannelIndex } }))
                    return true;
            }
            else
            {
                for (int i = 0; i < node->getProcessor()->getTotalNumInputChannels(); ++i)
                    if (i != inputChannelOfIndexToIgnore && graph.isConnected ({ output, { node->nodeId, i } }))
                        return true;
            }

            inputChannelOfIndexToIgnore = -1;
            ++stepIndexToSearchFrom;
        }

        return false;
    }

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (RenderSequenceBuilder)
};

//==============================================================================
AudioProcessorGraph::Connection::Connection (NodeAndChannel src, NodeAndChannel dst) noexcept
    : source (src), destination (dst)
{
}

bool AudioProcessorGraph::Connection::operator== (const Connection& other) const noexcept
{
    return source == other.source && destination == other.destination;
}

bool AudioProcessorGraph::Connection::operator!= (const Connection& c) const noexcept
{
    return ! operator== (c);
}

bool AudioProcessorGraph::Connection::operator< (const Connection& other) const noexcept
{
    if (source.nodeID != other.source.nodeID)
        return source.nodeID < other.source.nodeID;

    if (destination.nodeID != other.destination.nodeID)
        return destination.nodeID < other.destination.nodeID;

    if (source.channelIndex != other.source.channelIndex)
        return source.channelIndex < other.source.channelIndex;

    return destination.channelIndex < other.destination.channelIndex;
}

//==============================================================================
AudioProcessorGraph::Node::Node (NodeID nodeID, AudioProcessor* p) noexcept
    : nodeId (nodeID), processor (p)
{
    jassert (processor != nullptr);
}

void AudioProcessorGraph::Node::prepare (double newSampleRate, int newBlockSize,
                                         AudioProcessorGraph* graph, ProcessingPrecision precision)
{
    if (! isPrepared)
    {
        isPrepared = true;
        setParentGraph (graph);

        // try to align the precision of the processor and the graph
        processor->setProcessingPrecision (processor->supportsDoublePrecisionProcessing() ? precision
                                                                                          : singlePrecision);

        processor->setRateAndBufferSizeDetails (newSampleRate, newBlockSize);
        processor->prepareToPlay (newSampleRate, newBlockSize);
    }
}

void AudioProcessorGraph::Node::unprepare()
{
    if (isPrepared)
    {
        isPrepared = false;
        processor->releaseResources();
    }
}

void AudioProcessorGraph::Node::setParentGraph (AudioProcessorGraph* const graph) const
{
    if (auto* ioProc = dynamic_cast<AudioProcessorGraph::AudioGraphIOProcessor*> (processor.get()))
        ioProc->setParentGraph (graph);
}

//==============================================================================
struct AudioProcessorGraph::RenderSequenceFloat   : public GraphRenderSequence<float> {};
struct AudioProcessorGraph::RenderSequenceDouble  : public GraphRenderSequence<double> {};

//==============================================================================
AudioProcessorGraph::AudioProcessorGraph()
{
}

AudioProcessorGraph::~AudioProcessorGraph()
{
    clearRenderingSequence();
    clear();
}

const String AudioProcessorGraph::getName() const
{
    return "Audio Graph";
}

//==============================================================================
void AudioProcessorGraph::clear()
{
    nodes.clear();
    connections.clear();
    triggerAsyncUpdate();
}

AudioProcessorGraph::Node* AudioProcessorGraph::getNodeForId (NodeID nodeId) const
{
    for (auto* n : nodes)
        if (n->nodeId == nodeId)
            return n;

    return nullptr;
}

AudioProcessorGraph::Node* AudioProcessorGraph::addNode (AudioProcessor* newProcessor, NodeID nodeId)
{
    if (newProcessor == nullptr || newProcessor == this)
    {
        jassertfalse;
        return nullptr;
    }

    for (auto* n : nodes)
    {
        if (n->getProcessor() == newProcessor)
        {
            jassertfalse; // Cannot add the same object to the graph twice!
            return nullptr;
        }
    }

    if (nodeId == 0)
    {
        nodeId = ++lastNodeId;
    }
    else
    {
        // you can't add a node with an id that already exists in the graph..
        jassert (getNodeForId (nodeId) == nullptr);
        removeNode (nodeId);

        if (nodeId > lastNodeId)
            lastNodeId = nodeId;
    }

    newProcessor->setPlayHead (getPlayHead());

    auto* n = new Node (nodeId, newProcessor);
    nodes.add (n);

    if (isPrepared)
        triggerAsyncUpdate();

    n->setParentGraph (this);
    return n;
}

bool AudioProcessorGraph::removeNode (NodeID nodeId)
{
    disconnectNode (nodeId);

    for (int i = nodes.size(); --i >= 0;)
    {
        if (nodes.getUnchecked(i)->nodeId == nodeId)
        {
            nodes.remove (i);

            if (isPrepared)
                triggerAsyncUpdate();

            return true;
        }
    }

    return false;
}

bool AudioProcessorGraph::removeNode (Node* node)
{
    if (node != nullptr)
        return removeNode (node->nodeId);

    jassertfalse;
    return false;
}

//==============================================================================
bool AudioProcessorGraph::isConnected (const Connection& c) const noexcept
{
    sortConnections();

    size_t s = 0, e = connections.size();

    while (s < e)
    {
        if (c == connections[s])
            return true;

        auto halfway = (s + e) / 2;

        if (halfway == s)
            break;

        if (c < connections[halfway])
            e = halfway;
        else
            s = halfway;
    }

    return false;
}

bool AudioProcessorGraph::isConnected (NodeID srcID, NodeID destID) const noexcept
{
    for (auto&& c : connections)
        if (c.source.nodeID == srcID && c.destination.nodeID == destID)
            return true;

    return false;
}

bool AudioProcessorGraph::canConnect (const Connection& c) const
{
    if (c.source.channelIndex < 0
         || c.destination.channelIndex < 0
         || c.source.nodeID == c.destination.nodeID
         || c.destination.isMIDI() != c.source.isMIDI())
        return false;

    auto* source = getNodeForId (c.source.nodeID);

    if (source == nullptr
         || (! c.source.isMIDI() && c.source.channelIndex >= source->processor->getTotalNumOutputChannels())
         || (c.source.isMIDI() && ! source->processor->producesMidi()))
        return false;

    auto* dest = getNodeForId (c.destination.nodeID);

    if (dest == nullptr
         || (! c.destination.isMIDI() && c.destination.channelIndex >= dest->processor->getTotalNumInputChannels())
         || (c.destination.isMIDI() && ! dest->processor->acceptsMidi()))
        return false;

    return ! isConnected (c);
}

bool AudioProcessorGraph::addConnection (const Connection& c)
{
    if (! canConnect (c))
        return false;

    if (! (connectionsNeedSorting || connections.empty()))
        connectionsNeedSorting = (c < connections.back());

    connections.push_back (c);

    if (isPrepared)
        triggerAsyncUpdate();

    return true;
}

void AudioProcessorGraph::sortConnections() const noexcept
{
    if (connectionsNeedSorting)
    {
        connectionsNeedSorting = false;
        std::sort (connections.begin(), connections.end());
    }
}

bool AudioProcessorGraph::removeConnection (const Connection& c)
{
    auto pos = std::find (connections.begin(), connections.end(), c);

    if (pos != connections.end())
    {
        connections.erase (pos);

        if (isPrepared)
            triggerAsyncUpdate();

        return true;
    }

    return false;
}

bool AudioProcessorGraph::disconnectNode (NodeID nodeID)
{
    auto oldSize = connections.size();

    connections.erase (std::remove_if (connections.begin(), connections.end(),
                                      [nodeID] (const Connection& c) { return c.source.nodeID == nodeID || c.destination.nodeID == nodeID; }),
                       connections.end());

    return oldSize != connections.size();
}

bool AudioProcessorGraph::isConnectionLegal (const Connection& c) const
{
    if (auto* source = getNodeForId (c.source.nodeID))
        if (auto* dest = getNodeForId (c.destination.nodeID))
            return (c.source.isMIDI() ? source->processor->producesMidi()
                                      : isPositiveAndBelow (c.source.channelIndex, source->processor->getTotalNumOutputChannels()))
                && (c.destination.isMIDI() ? dest->processor->acceptsMidi()
                                           : isPositiveAndBelow (c.destination.channelIndex, dest->processor->getTotalNumInputChannels()));

    return false;
}

bool AudioProcessorGraph::removeIllegalConnections()
{
    auto oldSize = connections.size();

    connections.erase (std::remove_if (connections.begin(), connections.end(),
                                      [this] (const Connection& c) { return ! isConnectionLegal (c); }),
                       connections.end());

    return oldSize != connections.size();
}

//==============================================================================
void AudioProcessorGraph::clearRenderingSequence()
{
    ScopedPointer<RenderSequenceFloat> oldSequenceF;
    ScopedPointer<RenderSequenceDouble> oldSequenceD;

    {
        const ScopedLock sl (getCallbackLock());
        renderSequenceFloat.swapWith (oldSequenceF);
        renderSequenceDouble.swapWith (oldSequenceD);
    }
}

bool AudioProcessorGraph::isAnInputTo (NodeID src, NodeID dst, int recursionCheck) const
{
    if (recursionCheck > 0)
    {
        for (auto&& c : connections)
            if (c.destination.nodeID == dst
                 && (c.source.nodeID == src || isAnInputTo (src, c.source.nodeID, recursionCheck - 1)))
                return true;
    }

    return false;
}

bool AudioProcessorGraph::anyNodesNeedPreparing() const noexcept
{
    for (auto* node : nodes)
        if (! node->isPrepared)
            return true;

    return false;
}

void AudioProcessorGraph::buildRenderingSequence()
{
    ScopedPointer<RenderSequenceFloat>  newSequenceF (new RenderSequenceFloat());
    ScopedPointer<RenderSequenceDouble> newSequenceD (new RenderSequenceDouble());

    {
        MessageManagerLock mml;

        RenderSequenceBuilder<RenderSequenceFloat>  builderF (*this, *newSequenceF);
        RenderSequenceBuilder<RenderSequenceDouble> builderD (*this, *newSequenceD);
    }

    newSequenceF->prepareBuffers (getBlockSize());
    newSequenceD->prepareBuffers (getBlockSize());

    if (anyNodesNeedPreparing())
    {
        {
            const ScopedLock sl (getCallbackLock());
            renderSequenceFloat = nullptr;
            renderSequenceDouble = nullptr;
        }

        for (auto* node : nodes)
            node->prepare (getSampleRate(), getBlockSize(), this, getProcessingPrecision());
    }

    const ScopedLock sl (getCallbackLock());

    renderSequenceFloat.swapWith (newSequenceF);
    renderSequenceDouble.swapWith (newSequenceD);
}

void AudioProcessorGraph::handleAsyncUpdate()
{
    buildRenderingSequence();
}

//==============================================================================
void AudioProcessorGraph::prepareToPlay (double /*sampleRate*/, int estimatedSamplesPerBlock)
{
    if (renderSequenceFloat != nullptr)
        renderSequenceFloat->prepareBuffers (estimatedSamplesPerBlock);

    if (renderSequenceDouble != nullptr)
        renderSequenceDouble->prepareBuffers (estimatedSamplesPerBlock);

    clearRenderingSequence();
    buildRenderingSequence();

    isPrepared = true;
}

bool AudioProcessorGraph::supportsDoublePrecisionProcessing() const
{
    return true;
}

void AudioProcessorGraph::releaseResources()
{
    isPrepared = false;

    for (auto* n : nodes)
        n->unprepare();

    if (renderSequenceFloat != nullptr)
        renderSequenceFloat->releaseBuffers();

    if (renderSequenceDouble != nullptr)
        renderSequenceDouble->releaseBuffers();
}

void AudioProcessorGraph::reset()
{
    const ScopedLock sl (getCallbackLock());

    for (auto* n : nodes)
        n->getProcessor()->reset();
}

void AudioProcessorGraph::setNonRealtime (bool isProcessingNonRealtime) noexcept
{
    const ScopedLock sl (getCallbackLock());

    AudioProcessor::setNonRealtime (isProcessingNonRealtime);

    for (auto* n : nodes)
        n->getProcessor()->setNonRealtime (isProcessingNonRealtime);
}

void AudioProcessorGraph::setPlayHead (AudioPlayHead* audioPlayHead)
{
    const ScopedLock sl (getCallbackLock());

    AudioProcessor::setPlayHead (audioPlayHead);

    for (auto* n : nodes)
        n->getProcessor()->setPlayHead (audioPlayHead);
}

double AudioProcessorGraph::getTailLengthSeconds() const            { return 0; }
bool AudioProcessorGraph::acceptsMidi() const                       { return true; }
bool AudioProcessorGraph::producesMidi() const                      { return true; }
void AudioProcessorGraph::getStateInformation (juce::MemoryBlock&)  {}
void AudioProcessorGraph::setStateInformation (const void*, int)    {}

void AudioProcessorGraph::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    const ScopedLock sl (getCallbackLock());

    if (renderSequenceFloat != nullptr)
        renderSequenceFloat->perform (buffer, midiMessages);
}

void AudioProcessorGraph::processBlock (AudioBuffer<double>& buffer, MidiBuffer& midiMessages)
{
    const ScopedLock sl (getCallbackLock());

    if (renderSequenceDouble != nullptr)
        renderSequenceDouble->perform (buffer, midiMessages);
}

//==============================================================================
AudioProcessorGraph::AudioGraphIOProcessor::AudioGraphIOProcessor (const IODeviceType deviceType)
    : type (deviceType)
{
}

AudioProcessorGraph::AudioGraphIOProcessor::~AudioGraphIOProcessor()
{
}

const String AudioProcessorGraph::AudioGraphIOProcessor::getName() const
{
    switch (type)
    {
        case audioOutputNode:   return "Audio Output";
        case audioInputNode:    return "Audio Input";
        case midiOutputNode:    return "Midi Output";
        case midiInputNode:     return "Midi Input";
        default:                break;
    }

    return {};
}

void AudioProcessorGraph::AudioGraphIOProcessor::fillInPluginDescription (PluginDescription& d) const
{
    d.name = getName();
    d.uid = d.name.hashCode();
    d.category = "I/O devices";
    d.pluginFormatName = "Internal";
    d.manufacturerName = "ROLI Ltd.";
    d.version = "1.0";
    d.isInstrument = false;

    d.numInputChannels = getTotalNumInputChannels();

    if (type == audioOutputNode && graph != nullptr)
        d.numInputChannels = graph->getTotalNumInputChannels();

    d.numOutputChannels = getTotalNumOutputChannels();

    if (type == audioInputNode && graph != nullptr)
        d.numOutputChannels = graph->getTotalNumOutputChannels();
}

void AudioProcessorGraph::AudioGraphIOProcessor::prepareToPlay (double, int)
{
    jassert (graph != nullptr);
}

void AudioProcessorGraph::AudioGraphIOProcessor::releaseResources()
{
}

bool AudioProcessorGraph::AudioGraphIOProcessor::supportsDoublePrecisionProcessing() const
{
    return true;
}

template <typename FloatType, typename SequenceType>
static void processIOBlock (AudioProcessorGraph::AudioGraphIOProcessor& io, SequenceType& sequence,
                            AudioBuffer<FloatType>& buffer, MidiBuffer& midiMessages)
{
    switch (io.getType())
    {
        case AudioProcessorGraph::AudioGraphIOProcessor::audioOutputNode:
        {
            auto&& currentAudioOutputBuffer = sequence.currentAudioOutputBuffer;

            for (int i = jmin (currentAudioOutputBuffer.getNumChannels(), buffer.getNumChannels()); --i >= 0;)
                currentAudioOutputBuffer.addFrom (i, 0, buffer, i, 0, buffer.getNumSamples());

            break;
        }

        case AudioProcessorGraph::AudioGraphIOProcessor::audioInputNode:
        {
            auto* currentInputBuffer = sequence.currentAudioInputBuffer;

            for (int i = jmin (currentInputBuffer->getNumChannels(), buffer.getNumChannels()); --i >= 0;)
                buffer.copyFrom (i, 0, *currentInputBuffer, i, 0, buffer.getNumSamples());

            break;
        }

        case AudioProcessorGraph::AudioGraphIOProcessor::midiOutputNode:
            sequence.currentMidiOutputBuffer.addEvents (midiMessages, 0, buffer.getNumSamples(), 0);
            break;

        case AudioProcessorGraph::AudioGraphIOProcessor::midiInputNode:
            midiMessages.addEvents (*sequence.currentMidiInputBuffer, 0, buffer.getNumSamples(), 0);
            break;

        default:
            break;
    }
}

void AudioProcessorGraph::AudioGraphIOProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    jassert (graph != nullptr);
    processIOBlock (*this, *graph->renderSequenceFloat, buffer, midiMessages);
}

void AudioProcessorGraph::AudioGraphIOProcessor::processBlock (AudioBuffer<double>& buffer, MidiBuffer& midiMessages)
{
    jassert (graph != nullptr);
    processIOBlock (*this, *graph->renderSequenceDouble, buffer, midiMessages);
}

double AudioProcessorGraph::AudioGraphIOProcessor::getTailLengthSeconds() const
{
    return 0;
}

bool AudioProcessorGraph::AudioGraphIOProcessor::acceptsMidi() const
{
    return type == midiOutputNode;
}

bool AudioProcessorGraph::AudioGraphIOProcessor::producesMidi() const
{
    return type == midiInputNode;
}

bool AudioProcessorGraph::AudioGraphIOProcessor::isInput() const noexcept           { return type == audioInputNode  || type == midiInputNode; }
bool AudioProcessorGraph::AudioGraphIOProcessor::isOutput() const noexcept          { return type == audioOutputNode || type == midiOutputNode; }

bool AudioProcessorGraph::AudioGraphIOProcessor::hasEditor() const                  { return false; }
AudioProcessorEditor* AudioProcessorGraph::AudioGraphIOProcessor::createEditor()    { return nullptr; }

int AudioProcessorGraph::AudioGraphIOProcessor::getNumPrograms()                    { return 0; }
int AudioProcessorGraph::AudioGraphIOProcessor::getCurrentProgram()                 { return 0; }
void AudioProcessorGraph::AudioGraphIOProcessor::setCurrentProgram (int)            { }

const String AudioProcessorGraph::AudioGraphIOProcessor::getProgramName (int)       { return {}; }
void AudioProcessorGraph::AudioGraphIOProcessor::changeProgramName (int, const String&) {}

void AudioProcessorGraph::AudioGraphIOProcessor::getStateInformation (juce::MemoryBlock&) {}
void AudioProcessorGraph::AudioGraphIOProcessor::setStateInformation (const void*, int) {}

void AudioProcessorGraph::AudioGraphIOProcessor::setParentGraph (AudioProcessorGraph* const newGraph)
{
    graph = newGraph;

    if (graph != nullptr)
    {
        setPlayConfigDetails (type == audioOutputNode ? graph->getTotalNumOutputChannels() : 0,
                              type == audioInputNode  ? graph->getTotalNumInputChannels()  : 0,
                              getSampleRate(),
                              getBlockSize());

        updateHostDisplay();
    }
}
