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

const int AudioProcessorGraph::midiChannelIndex = 0x1000;


//==============================================================================
struct GraphConnectionSorter
{
    static int compareElements (const AudioProcessorGraph::Connection* first,
                                const AudioProcessorGraph::Connection* second) noexcept
    {
        if (first->sourceNodeId < second->sourceNodeId)                return -1;
        if (first->sourceNodeId > second->sourceNodeId)                return 1;
        if (first->destNodeId < second->destNodeId)                    return -1;
        if (first->destNodeId > second->destNodeId)                    return 1;
        if (first->sourceChannelIndex < second->sourceChannelIndex)    return -1;
        if (first->sourceChannelIndex > second->sourceChannelIndex)    return 1;
        if (first->destChannelIndex < second->destChannelIndex)        return -1;
        if (first->destChannelIndex > second->destChannelIndex)        return 1;

        return 0;
    }
};

//==============================================================================
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
            MidiBuffer endMIDI;
            endMIDI.addEvents (midiMessages, maxSamples, numSamples, -maxSamples);

            {
                AudioBuffer<FloatType> startAudio (buffer.getArrayOfWritePointers(), buffer.getNumChannels(), maxSamples);
                midiMessages.clear (maxSamples, numSamples);
                perform (startAudio, midiMessages);
            }

            AudioBuffer<FloatType> endAudio (buffer.getArrayOfWritePointers(), buffer.getNumChannels(), maxSamples, numSamples - maxSamples);
            perform (endAudio, endMIDI);
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
            {
                buffer.clear();
            }
            else
            {
                ScopedLock lock (processor.getCallbackLock());

                callProcess (buffer, c.midiBuffers[midiBufferToUse]);
            }
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

        nodeIds.add (zeroNodeID()); // first buffer is read-only zeros
        channels.add (0);

        midiNodeIds.add (zeroNodeID());

        for (int i = 0; i < orderedNodes.size(); ++i)
        {
            createRenderingOpsForNode (*orderedNodes.getUnchecked(i), i);
            markAnyUnusedBuffersAsFree (i);
        }

        graph.setLatencySamples (totalLatency);

        s.numBuffersNeeded = nodeIds.size();
        s.numMidiBuffersNeeded = midiNodeIds.size();
    }

    //==============================================================================
    typedef AudioProcessorGraph::NodeID NodeID;

    AudioProcessorGraph& graph;
    RenderSequence& sequence;

    Array<AudioProcessorGraph::Node*> orderedNodes;
    Array<int> channels;
    Array<NodeID> nodeIds, midiNodeIds;

    static NodeID zeroNodeID()       { return 0xfffffffe; }
    static NodeID freeNodeID()       { return 0xffffffff; }
    static NodeID anonymousNodeID()  { return 0xfffffffd; }

    static bool isNodeBusy (NodeID nodeID) noexcept     { return nodeID != freeNodeID() && nodeID != zeroNodeID(); }

    Array<NodeID> nodeDelayIDs;
    Array<int> nodeDelays;
    int totalLatency = 0;

    int getNodeDelay (NodeID nodeID) const        { return nodeDelays [nodeDelayIDs.indexOf (nodeID)]; }

    void setNodeDelay (NodeID nodeID, const int latency)
    {
        const int index = nodeDelayIDs.indexOf (nodeID);

        if (index >= 0)
        {
            nodeDelays.set (index, latency);
        }
        else
        {
            nodeDelayIDs.add (nodeID);
            nodeDelays.add (latency);
        }
    }

    int getInputLatencyForNode (NodeID nodeID) const
    {
        int maxLatency = 0;

        for (auto* c : graph.getConnections())
            if (c->destNodeId == nodeID)
                maxLatency = jmax (maxLatency, getNodeDelay (c->sourceNodeId));

        return maxLatency;
    }

    //==============================================================================
    // Fast lookup table for checking which nodes are inputs to others.
    struct NodeInputLookupTable
    {
        template <typename ConnectionList>
        explicit NodeInputLookupTable (const ConnectionList& connections)
        {
            for (auto* c : connections)
                entries[c->destNodeId].add (c->sourceNodeId);
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

    void createRenderingOpsForNode (AudioProcessorGraph::Node& node, const int ourRenderingIndex)
    {
        auto& processor = *node.getProcessor();
        auto numIns  = processor.getTotalNumInputChannels();
        auto numOuts = processor.getTotalNumOutputChannels();
        auto totalChans = jmax (numIns, numOuts);

        Array<int> audioChannelsToUse;
        int midiBufferToUse = -1;

        int maxLatency = getInputLatencyForNode (node.nodeId);

        for (int inputChan = 0; inputChan < numIns; ++inputChan)
        {
            // get a list of all the inputs to this node
            Array<NodeID> sourceNodes;
            Array<int> sourceOutputChans;

            for (auto* c : graph.getConnections())
            {
                if (c->destNodeId == node.nodeId && c->destChannelIndex == inputChan)
                {
                    sourceNodes.add (c->sourceNodeId);
                    sourceOutputChans.add (c->sourceChannelIndex);
                }
            }

            int bufIndex = -1;

            if (sourceNodes.isEmpty())
            {
                // unconnected input channel

                if (inputChan >= numOuts)
                {
                    bufIndex = getReadOnlyEmptyBuffer();
                    jassert (bufIndex >= 0);
                }
                else
                {
                    bufIndex = getFreeBuffer (false);
                    sequence.addClearChannelOp (bufIndex);
                }
            }
            else if (sourceNodes.size() == 1)
            {
                // channel with a straightforward single input..
                auto srcNode = sourceNodes.getUnchecked(0);
                auto srcChan = sourceOutputChans.getUnchecked(0);

                bufIndex = getBufferContaining (srcNode, srcChan);

                if (bufIndex < 0)
                {
                    // if not found, this is probably a feedback loop
                    bufIndex = getReadOnlyEmptyBuffer();
                    jassert (bufIndex >= 0);
                }

                if (inputChan < numOuts
                     && isBufferNeededLater (ourRenderingIndex,
                                             inputChan,
                                             srcNode, srcChan))
                {
                    // can't mess up this channel because it's needed later by another node, so we
                    // need to use a copy of it..
                    const int newFreeBuffer = getFreeBuffer (false);

                    sequence.addCopyChannelOp (bufIndex, newFreeBuffer);

                    bufIndex = newFreeBuffer;
                }

                const int nodeDelay = getNodeDelay (srcNode);

                if (nodeDelay < maxLatency)
                    sequence.addDelayChannelOp (bufIndex, maxLatency - nodeDelay);
            }
            else
            {
                // channel with a mix of several inputs..

                // try to find a re-usable channel from our inputs..
                int reusableInputIndex = -1;

                for (int i = 0; i < sourceNodes.size(); ++i)
                {
                    const int sourceBufIndex = getBufferContaining (sourceNodes.getUnchecked(i),
                                                                    sourceOutputChans.getUnchecked(i));

                    if (sourceBufIndex >= 0
                        && ! isBufferNeededLater (ourRenderingIndex,
                                                  inputChan,
                                                  sourceNodes.getUnchecked(i),
                                                  sourceOutputChans.getUnchecked(i)))
                    {
                        // we've found one of our input chans that can be re-used..
                        reusableInputIndex = i;
                        bufIndex = sourceBufIndex;

                        const int nodeDelay = getNodeDelay (sourceNodes.getUnchecked (i));
                        if (nodeDelay < maxLatency)
                            sequence.addDelayChannelOp (sourceBufIndex, maxLatency - nodeDelay);

                        break;
                    }
                }

                if (reusableInputIndex < 0)
                {
                    // can't re-use any of our input chans, so get a new one and copy everything into it..
                    bufIndex = getFreeBuffer (false);
                    jassert (bufIndex != 0);

                    markBufferAsContaining (bufIndex, anonymousNodeID(), 0);

                    const int srcIndex = getBufferContaining (sourceNodes.getUnchecked (0),
                                                              sourceOutputChans.getUnchecked (0));
                    if (srcIndex < 0)
                        sequence.addClearChannelOp (bufIndex);  // if not found, this is probably a feedback loop
                    else
                        sequence.addCopyChannelOp (srcIndex, bufIndex);

                    reusableInputIndex = 0;
                    const int nodeDelay = getNodeDelay (sourceNodes.getFirst());

                    if (nodeDelay < maxLatency)
                        sequence.addDelayChannelOp (bufIndex, maxLatency - nodeDelay);
                }

                for (int j = 0; j < sourceNodes.size(); ++j)
                {
                    if (j != reusableInputIndex)
                    {
                        int srcIndex = getBufferContaining (sourceNodes.getUnchecked(j),
                                                            sourceOutputChans.getUnchecked(j));
                        if (srcIndex >= 0)
                        {
                            const int nodeDelay = getNodeDelay (sourceNodes.getUnchecked (j));

                            if (nodeDelay < maxLatency)
                            {
                                if (! isBufferNeededLater (ourRenderingIndex, inputChan,
                                                           sourceNodes.getUnchecked(j),
                                                           sourceOutputChans.getUnchecked(j)))
                                {
                                    sequence.addDelayChannelOp (srcIndex, maxLatency - nodeDelay);
                                }
                                else // buffer is reused elsewhere, can't be delayed
                                {
                                    const int bufferToDelay = getFreeBuffer (false);
                                    sequence.addCopyChannelOp (srcIndex, bufferToDelay);
                                    sequence.addDelayChannelOp (bufferToDelay, maxLatency - nodeDelay);
                                    srcIndex = bufferToDelay;
                                }
                            }

                            sequence.addAddChannelOp (srcIndex, bufIndex);
                        }
                    }
                }
            }

            jassert (bufIndex >= 0);
            audioChannelsToUse.add (bufIndex);

            if (inputChan < numOuts)
                markBufferAsContaining (bufIndex, node.nodeId, inputChan);
        }

        for (int outputChan = numIns; outputChan < numOuts; ++outputChan)
        {
            const int bufIndex = getFreeBuffer (false);
            jassert (bufIndex != 0);
            audioChannelsToUse.add (bufIndex);

            markBufferAsContaining (bufIndex, node.nodeId, outputChan);
        }

        // Now the same thing for midi..
        Array<NodeID> midiSourceNodes;

        for (int i = graph.getNumConnections(); --i >= 0;)
        {
            auto* c = graph.getConnection (i);

            if (c->destNodeId == node.nodeId && c->destChannelIndex == AudioProcessorGraph::midiChannelIndex)
                midiSourceNodes.add (c->sourceNodeId);
        }

        if (midiSourceNodes.isEmpty())
        {
            // No midi inputs..
            midiBufferToUse = getFreeBuffer (true); // need to pick a buffer even if the processor doesn't use midi

            if (processor.acceptsMidi() || processor.producesMidi())
                sequence.addClearMidiBufferOp (midiBufferToUse);
        }
        else if (midiSourceNodes.size() == 1)
        {
            // One midi input..
            midiBufferToUse = getBufferContaining (midiSourceNodes.getUnchecked(0),
                                                   AudioProcessorGraph::midiChannelIndex);

            if (midiBufferToUse >= 0)
            {
                if (isBufferNeededLater (ourRenderingIndex,
                                         AudioProcessorGraph::midiChannelIndex,
                                         midiSourceNodes.getUnchecked(0),
                                         AudioProcessorGraph::midiChannelIndex))
                {
                    // can't mess up this channel because it's needed later by another node, so we
                    // need to use a copy of it..
                    const int newFreeBuffer = getFreeBuffer (true);
                    sequence.addCopyMidiBufferOp (midiBufferToUse, newFreeBuffer);
                    midiBufferToUse = newFreeBuffer;
                }
            }
            else
            {
                // probably a feedback loop, so just use an empty one..
                midiBufferToUse = getFreeBuffer (true); // need to pick a buffer even if the processor doesn't use midi
            }
        }
        else
        {
            // More than one midi input being mixed..
            int reusableInputIndex = -1;

            for (int i = 0; i < midiSourceNodes.size(); ++i)
            {
                const int sourceBufIndex = getBufferContaining (midiSourceNodes.getUnchecked(i),
                                                                AudioProcessorGraph::midiChannelIndex);

                if (sourceBufIndex >= 0
                     && ! isBufferNeededLater (ourRenderingIndex,
                                               AudioProcessorGraph::midiChannelIndex,
                                               midiSourceNodes.getUnchecked(i),
                                               AudioProcessorGraph::midiChannelIndex))
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
                midiBufferToUse = getFreeBuffer (true);
                jassert (midiBufferToUse >= 0);

                const int srcIndex = getBufferContaining (midiSourceNodes.getUnchecked(0),
                                                          AudioProcessorGraph::midiChannelIndex);
                if (srcIndex >= 0)
                    sequence.addCopyMidiBufferOp (srcIndex, midiBufferToUse);
                else
                    sequence.addClearMidiBufferOp (midiBufferToUse);

                reusableInputIndex = 0;
            }

            for (int j = 0; j < midiSourceNodes.size(); ++j)
            {
                if (j != reusableInputIndex)
                {
                    const int srcIndex = getBufferContaining (midiSourceNodes.getUnchecked(j),
                                                              AudioProcessorGraph::midiChannelIndex);
                    if (srcIndex >= 0)
                        sequence.addAddMidiBufferOp (srcIndex, midiBufferToUse);
                }
            }
        }

        if (processor.producesMidi())
            markBufferAsContaining (midiBufferToUse, node.nodeId,
                                    AudioProcessorGraph::midiChannelIndex);

        setNodeDelay (node.nodeId, maxLatency + processor.getLatencySamples());

        if (numOuts == 0)
            totalLatency = maxLatency;

        sequence.addProcessOp (&node, audioChannelsToUse, totalChans, midiBufferToUse);
    }

    //==============================================================================
    int getFreeBuffer (const bool forMidi)
    {
        if (forMidi)
        {
            for (int i = 1; i < midiNodeIds.size(); ++i)
                if (midiNodeIds.getUnchecked(i) == freeNodeID())
                    return i;

            midiNodeIds.add (freeNodeID());
            return midiNodeIds.size() - 1;
        }

        for (int i = 1; i < nodeIds.size(); ++i)
            if (nodeIds.getUnchecked(i) == freeNodeID())
                return i;

        nodeIds.add (freeNodeID());
        channels.add (0);
        return nodeIds.size() - 1;
    }

    int getReadOnlyEmptyBuffer() const noexcept
    {
        return 0;
    }

    int getBufferContaining (NodeID nodeId, int outputChannel) const noexcept
    {
        if (outputChannel == AudioProcessorGraph::midiChannelIndex)
        {
            for (int i = midiNodeIds.size(); --i >= 0;)
                if (midiNodeIds.getUnchecked(i) == nodeId)
                    return i;
        }
        else
        {
            for (int i = nodeIds.size(); --i >= 0;)
                if (nodeIds.getUnchecked(i) == nodeId
                     && channels.getUnchecked(i) == outputChannel)
                    return i;
        }

        return -1;
    }

    void markAnyUnusedBuffersAsFree (const int stepIndex)
    {
        for (int i = 0; i < nodeIds.size(); ++i)
        {
            if (isNodeBusy (nodeIds.getUnchecked(i))
                 && ! isBufferNeededLater (stepIndex, -1,
                                           nodeIds.getUnchecked(i),
                                           channels.getUnchecked(i)))
            {
                nodeIds.set (i, freeNodeID());
            }
        }

        for (int i = 0; i < midiNodeIds.size(); ++i)
        {
            if (isNodeBusy (midiNodeIds.getUnchecked(i))
                 && ! isBufferNeededLater (stepIndex, -1,
                                           midiNodeIds.getUnchecked(i),
                                           AudioProcessorGraph::midiChannelIndex))
            {
                midiNodeIds.set (i, freeNodeID());
            }
        }
    }

    bool isBufferNeededLater (int stepIndexToSearchFrom,
                              int inputChannelOfIndexToIgnore,
                              NodeID nodeId,
                              int outputChanIndex) const
    {
        while (stepIndexToSearchFrom < orderedNodes.size())
        {
            auto* node = orderedNodes.getUnchecked (stepIndexToSearchFrom);

            if (outputChanIndex == AudioProcessorGraph::midiChannelIndex)
            {
                if (inputChannelOfIndexToIgnore != AudioProcessorGraph::midiChannelIndex
                     && graph.getConnectionBetween (nodeId, AudioProcessorGraph::midiChannelIndex,
                                                    node->nodeId, AudioProcessorGraph::midiChannelIndex) != nullptr)
                    return true;
            }
            else
            {
                for (int i = 0; i < node->getProcessor()->getTotalNumInputChannels(); ++i)
                    if (i != inputChannelOfIndexToIgnore
                         && graph.getConnectionBetween (nodeId, outputChanIndex,
                                                        node->nodeId, i) != nullptr)
                        return true;
            }

            inputChannelOfIndexToIgnore = -1;
            ++stepIndexToSearchFrom;
        }

        return false;
    }

    void markBufferAsContaining (int bufferNum, NodeID nodeId, int outputIndex)
    {
        if (outputIndex == AudioProcessorGraph::midiChannelIndex)
        {
            jassert (bufferNum > 0 && bufferNum < midiNodeIds.size());

            midiNodeIds.set (bufferNum, nodeId);
        }
        else
        {
            jassert (bufferNum >= 0 && bufferNum < nodeIds.size());

            nodeIds.set (bufferNum, nodeId);
            channels.set (bufferNum, outputIndex);
        }
    }

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (RenderSequenceBuilder)
};

//==============================================================================
AudioProcessorGraph::Connection::Connection (NodeID sourceID, int sourceChannel,
                                             NodeID destID, int destChannel) noexcept
    : sourceNodeId (sourceID), sourceChannelIndex (sourceChannel),
      destNodeId (destID), destChannelIndex (destChannel)
{
}

//==============================================================================
AudioProcessorGraph::Node::Node (NodeID nodeID, AudioProcessor* p) noexcept
    : nodeId (nodeID), processor (p)
{
    jassert (processor != nullptr);
}

void AudioProcessorGraph::Node::prepare (const double newSampleRate, const int newBlockSize,
                                         AudioProcessorGraph* const graph, ProcessingPrecision precision)
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
const AudioProcessorGraph::Connection* AudioProcessorGraph::getConnectionBetween (NodeID sourceNodeId, int sourceChannelIndex,
                                                                                  NodeID destNodeId, int destChannelIndex) const
{
    const Connection c (sourceNodeId, sourceChannelIndex, destNodeId, destChannelIndex);
    GraphConnectionSorter sorter;
    return connections [connections.indexOfSorted (sorter, &c)];
}

bool AudioProcessorGraph::isConnected (NodeID possibleSourceNodeId, NodeID possibleDestNodeId) const
{
    for (auto* c : connections)
        if (c->sourceNodeId == possibleSourceNodeId && c->destNodeId == possibleDestNodeId)
            return true;

    return false;
}

bool AudioProcessorGraph::canConnect (NodeID sourceNodeId, int sourceChannelIndex,
                                      NodeID destNodeId, int destChannelIndex) const
{
    if (sourceChannelIndex < 0
         || destChannelIndex < 0
         || sourceNodeId == destNodeId
         || (destChannelIndex == midiChannelIndex) != (sourceChannelIndex == midiChannelIndex))
        return false;

    auto* source = getNodeForId (sourceNodeId);

    if (source == nullptr
         || (sourceChannelIndex != midiChannelIndex && sourceChannelIndex >= source->processor->getTotalNumOutputChannels())
         || (sourceChannelIndex == midiChannelIndex && ! source->processor->producesMidi()))
        return false;

    auto* dest = getNodeForId (destNodeId);

    if (dest == nullptr
         || (destChannelIndex != midiChannelIndex && destChannelIndex >= dest->processor->getTotalNumInputChannels())
         || (destChannelIndex == midiChannelIndex && ! dest->processor->acceptsMidi()))
        return false;

    return getConnectionBetween (sourceNodeId, sourceChannelIndex,
                                 destNodeId, destChannelIndex) == nullptr;
}

bool AudioProcessorGraph::addConnection (NodeID sourceNodeId, int sourceChannelIndex,
                                         NodeID destNodeId, int destChannelIndex)
{
    if (! canConnect (sourceNodeId, sourceChannelIndex, destNodeId, destChannelIndex))
        return false;

    GraphConnectionSorter sorter;
    connections.addSorted (sorter, new Connection (sourceNodeId, sourceChannelIndex,
                                                   destNodeId, destChannelIndex));

    if (isPrepared)
        triggerAsyncUpdate();

    return true;
}

void AudioProcessorGraph::removeConnection (const int index)
{
    connections.remove (index);

    if (isPrepared)
        triggerAsyncUpdate();
}

bool AudioProcessorGraph::removeConnection (NodeID sourceNodeId, int sourceChannelIndex,
                                            NodeID destNodeId, int destChannelIndex)
{
    bool doneAnything = false;

    for (int i = connections.size(); --i >= 0;)
    {
        auto* c = connections.getUnchecked(i);

        if (c->sourceNodeId == sourceNodeId
             && c->destNodeId == destNodeId
             && c->sourceChannelIndex == sourceChannelIndex
             && c->destChannelIndex == destChannelIndex)
        {
            removeConnection (i);
            doneAnything = true;
        }
    }

    return doneAnything;
}

bool AudioProcessorGraph::disconnectNode (NodeID nodeId)
{
    bool doneAnything = false;

    for (int i = connections.size(); --i >= 0;)
    {
        auto* c = connections.getUnchecked(i);

        if (c->sourceNodeId == nodeId || c->destNodeId == nodeId)
        {
            removeConnection (i);
            doneAnything = true;
        }
    }

    return doneAnything;
}

bool AudioProcessorGraph::isConnectionLegal (const Connection* const c) const
{
    jassert (c != nullptr);

    if (auto* source = getNodeForId (c->sourceNodeId))
        if (auto* dest = getNodeForId (c->destNodeId))
            return (c->sourceChannelIndex != midiChannelIndex ? isPositiveAndBelow (c->sourceChannelIndex, source->processor->getTotalNumOutputChannels())
                                                              : source->processor->producesMidi())
                && (c->destChannelIndex   != midiChannelIndex ? isPositiveAndBelow (c->destChannelIndex, dest->processor->getTotalNumInputChannels())
                                                              : dest->processor->acceptsMidi());

    return false;
}

bool AudioProcessorGraph::removeIllegalConnections()
{
    bool doneAnything = false;

    for (int i = connections.size(); --i >= 0;)
    {
        if (! isConnectionLegal (connections.getUnchecked(i)))
        {
            removeConnection (i);
            doneAnything = true;
        }
    }

    return doneAnything;
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
        for (auto* c : connections)
            if (c->destNodeId == dst
                 && (c->sourceNodeId == src || isAnInputTo (src, c->sourceNodeId, recursionCheck - 1)))
                return true;
    }

    return false;
}

void AudioProcessorGraph::buildRenderingSequence()
{
    ScopedPointer<RenderSequenceFloat> newSequenceF;
    ScopedPointer<RenderSequenceDouble> newSequenceD;

    {
        MessageManagerLock mml;

        for (auto* node : nodes)
            node->prepare (getSampleRate(), getBlockSize(), this, getProcessingPrecision());

        newSequenceF = new RenderSequenceFloat();
        newSequenceD = new RenderSequenceDouble();

        RenderSequenceBuilder<RenderSequenceFloat> builderF (*this, *newSequenceF);
        RenderSequenceBuilder<RenderSequenceDouble> builderD (*this, *newSequenceD);

        newSequenceF->prepareBuffers (getBlockSize());
        newSequenceD->prepareBuffers (getBlockSize());
    }

    {
        // swap over to the new rendering sequence..
        const ScopedLock sl (getCallbackLock());
        renderSequenceFloat.swapWith (newSequenceF);
        renderSequenceDouble.swapWith (newSequenceD);
    }
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
    if (renderSequenceFloat != nullptr)
        renderSequenceFloat->perform (buffer, midiMessages);
}

void AudioProcessorGraph::processBlock (AudioBuffer<double>& buffer, MidiBuffer& midiMessages)
{
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
