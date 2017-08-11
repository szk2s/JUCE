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
template <typename FloatType, typename Impl> struct FloatDoubleUtil {};
template <typename Tag, typename Type> struct FloatDoubleType {};

template <typename Tag>
struct FloatAndDoubleComposition
{
    typedef typename FloatDoubleType<Tag, float>::Type  FloatType;
    typedef typename FloatDoubleType<Tag, double>::Type DoubleType;

    template <typename FloatingType>
    inline typename FloatDoubleType<Tag, FloatingType>::Type& get() noexcept
    {
        return FloatDoubleUtil<FloatingType, FloatAndDoubleComposition<Tag>>::get (*this);
    }

    FloatType floatVersion;
    DoubleType doubleVersion;
};

template <typename Impl> struct FloatDoubleUtil<float,  Impl> { static inline typename Impl::FloatType&  get (Impl& i) noexcept { return i.floatVersion; } };
template <typename Impl> struct FloatDoubleUtil<double, Impl> { static inline typename Impl::DoubleType& get (Impl& i) noexcept { return i.doubleVersion; } };

struct FloatPlaceholder;

template <typename FloatingType> struct FloatDoubleType<HeapBlock<FloatPlaceholder>,    FloatingType>  { typedef HeapBlock<FloatingType> Type; };
template <typename FloatingType> struct FloatDoubleType<HeapBlock<FloatPlaceholder*>,   FloatingType>  { typedef HeapBlock<FloatingType*> Type; };
template <typename FloatingType> struct FloatDoubleType<AudioBuffer<FloatPlaceholder>,  FloatingType>  { typedef AudioBuffer<FloatingType> Type; };
template <typename FloatingType> struct FloatDoubleType<AudioBuffer<FloatPlaceholder>*, FloatingType>  { typedef AudioBuffer<FloatingType>* Type; };


//==============================================================================
struct ConnectionSorter
{
    static int compareElements (const AudioProcessorGraph::Connection* const first,
                                const AudioProcessorGraph::Connection* const second) noexcept
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
struct AudioProcessorGraph::RenderSequence
{
    RenderSequence (AudioProcessorGraph& g)
    {
    }

    template <typename FloatType>
    void perform (AudioBuffer<FloatType>& buffer, const OwnedArray<MidiBuffer>& sharedMidiBuffers, int numSamples)
    {
        for (auto* op : renderOps)
            op->perform (buffer, sharedMidiBuffers, numSamples);
    }

    int numBuffersNeeded = 0, numMidiBuffersNeeded = 0;

    void addClearChannelOp (int channel)                    { createOp<ClearChannelOp> (channel); }
    void addCopyChannelOp (int srcIndex, int dstIndex)      { createOp<CopyChannelOp> (srcIndex, dstIndex); }
    void addAddChannelOp (int srcIndex, int dstIndex)       { createOp<AddChannelOp> (srcIndex, dstIndex); }
    void addClearMidiBufferOp (int index)                   { createOp<ClearMidiBufferOp> (index); }
    void addCopyMidiBufferOp (int srcIndex, int dstIndex)   { createOp<CopyMidiBufferOp> (srcIndex, dstIndex); }
    void addAddMidiBufferOp (int srcIndex, int dstIndex)    { createOp<AddMidiBufferOp> (srcIndex, dstIndex); }
    void addDelayChannelOp (int chan, int delaySize)        { createOp<DelayChannelOp> (chan, delaySize); }

    void addProcessBufferOp (const AudioProcessorGraph::Node::Ptr& node,
                             const Array<int>& audioChannelsUsed, int totalNumChans, int midiBuffer)
    {
        createOp<ProcessBufferOp> (node, audioChannelsUsed, totalNumChans, midiBuffer);
    }

private:
    struct RenderingOpBase
    {
        RenderingOpBase() noexcept {}
        virtual ~RenderingOpBase() {}

        virtual void perform (AudioBuffer<float>&,  const OwnedArray<MidiBuffer>&, int numSamples) = 0;
        virtual void perform (AudioBuffer<double>&, const OwnedArray<MidiBuffer>&, int numSamples) = 0;

        JUCE_LEAK_DETECTOR (RenderingOpBase)
    };

    OwnedArray<RenderingOpBase> renderOps;

    template <typename OpType, typename... Args>
    void createOp (Args... args) { renderOps.add (new OpType (args...)); }

    // use CRTP
    template <class Child>
    struct RenderingOp  : public RenderingOpBase
    {
        void perform (AudioBuffer<float>& sharedBufferChans,
                      const OwnedArray<MidiBuffer>& sharedMidiBuffers,
                      const int numSamples) override
        {
            static_cast<Child*> (this)->perform (sharedBufferChans, sharedMidiBuffers, numSamples);
        }

        void perform (AudioBuffer<double>& sharedBufferChans,
                      const OwnedArray<MidiBuffer>& sharedMidiBuffers,
                      const int numSamples) override
        {
            static_cast<Child*> (this)->perform (sharedBufferChans, sharedMidiBuffers, numSamples);
        }
    };

    //==============================================================================
    struct ClearChannelOp  : public RenderingOp<ClearChannelOp>
    {
        ClearChannelOp (const int channel) noexcept  : channelNum (channel)  {}

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>& sharedBufferChans, const OwnedArray<MidiBuffer>&, const int numSamples)
        {
            sharedBufferChans.clear (channelNum, 0, numSamples);
        }

        const int channelNum;

        JUCE_DECLARE_NON_COPYABLE (ClearChannelOp)
    };

    //==============================================================================
    struct CopyChannelOp  : public RenderingOp<CopyChannelOp>
    {
        CopyChannelOp (const int srcChan, const int dstChan) noexcept
            : srcChannelNum (srcChan), dstChannelNum (dstChan)
        {}

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>& sharedBufferChans, const OwnedArray<MidiBuffer>&, const int numSamples)
        {
            sharedBufferChans.copyFrom (dstChannelNum, 0, sharedBufferChans, srcChannelNum, 0, numSamples);
        }

        const int srcChannelNum, dstChannelNum;

        JUCE_DECLARE_NON_COPYABLE (CopyChannelOp)
    };

    //==============================================================================
    struct AddChannelOp  : public RenderingOp<AddChannelOp>
    {
        AddChannelOp (const int srcChan, const int dstChan) noexcept
            : srcChannelNum (srcChan), dstChannelNum (dstChan)
        {}

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>& sharedBufferChans, const OwnedArray<MidiBuffer>&, const int numSamples)
        {
            sharedBufferChans.addFrom (dstChannelNum, 0, sharedBufferChans, srcChannelNum, 0, numSamples);
        }

        const int srcChannelNum, dstChannelNum;

        JUCE_DECLARE_NON_COPYABLE (AddChannelOp)
    };

    //==============================================================================
    struct ClearMidiBufferOp  : public RenderingOp<ClearMidiBufferOp>
    {
        ClearMidiBufferOp (const int buffer) noexcept  : bufferNum (buffer)  {}

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>&, const OwnedArray<MidiBuffer>& sharedMidiBuffers, const int)
        {
            sharedMidiBuffers.getUnchecked (bufferNum)->clear();
        }

        const int bufferNum;

        JUCE_DECLARE_NON_COPYABLE (ClearMidiBufferOp)
    };

    //==============================================================================
    struct CopyMidiBufferOp  : public RenderingOp<CopyMidiBufferOp>
    {
        CopyMidiBufferOp (const int srcBuffer, const int dstBuffer) noexcept
            : srcBufferNum (srcBuffer), dstBufferNum (dstBuffer)
        {}

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>&, const OwnedArray<MidiBuffer>& sharedMidiBuffers, const int)
        {
            *sharedMidiBuffers.getUnchecked (dstBufferNum) = *sharedMidiBuffers.getUnchecked (srcBufferNum);
        }

        const int srcBufferNum, dstBufferNum;

        JUCE_DECLARE_NON_COPYABLE (CopyMidiBufferOp)
    };

    //==============================================================================
    struct AddMidiBufferOp  : public RenderingOp<AddMidiBufferOp>
    {
        AddMidiBufferOp (const int srcBuffer, const int dstBuffer)
            : srcBufferNum (srcBuffer), dstBufferNum (dstBuffer)
        {}

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>&, const OwnedArray<MidiBuffer>& sharedMidiBuffers, const int numSamples)
        {
            sharedMidiBuffers.getUnchecked (dstBufferNum)
                ->addEvents (*sharedMidiBuffers.getUnchecked (srcBufferNum), 0, numSamples, 0);
        }

        const int srcBufferNum, dstBufferNum;

        JUCE_DECLARE_NON_COPYABLE (AddMidiBufferOp)
    };

    //==============================================================================
    struct DelayChannelOp  : public RenderingOp<DelayChannelOp>
    {
        DelayChannelOp (const int chan, const int delaySize)
            : channel (chan),
              bufferSize (delaySize + 1),
              readIndex (0), writeIndex (delaySize)
        {
            buffer.floatVersion. calloc ((size_t) bufferSize);
            buffer.doubleVersion.calloc ((size_t) bufferSize);
        }

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>& sharedBufferChans, const OwnedArray<MidiBuffer>&, const int numSamples)
        {
            auto* data = sharedBufferChans.getWritePointer (channel, 0);
            auto& block = buffer.get<FloatType>();

            for (int i = numSamples; --i >= 0;)
            {
                block [writeIndex] = *data;
                *data++ = block [readIndex];

                if (++readIndex  >= bufferSize) readIndex = 0;
                if (++writeIndex >= bufferSize) writeIndex = 0;
            }
        }

    private:
        FloatAndDoubleComposition<HeapBlock<FloatPlaceholder>> buffer;
        const int channel, bufferSize;
        int readIndex, writeIndex;

        JUCE_DECLARE_NON_COPYABLE (DelayChannelOp)
    };

    //==============================================================================
    struct ProcessBufferOp   : public RenderingOp<ProcessBufferOp>
    {
        ProcessBufferOp (const AudioProcessorGraph::Node::Ptr& n,
                         const Array<int>& audioChannelsUsed,
                         const int totalNumChans,
                         const int midiBuffer)
            : node (n),
              processor (n->getProcessor()),
              audioChannelsToUse (audioChannelsUsed),
              totalChans (jmax (1, totalNumChans)),
              midiBufferToUse (midiBuffer)
        {
            audioChannels.floatVersion. calloc ((size_t) totalChans);
            audioChannels.doubleVersion.calloc ((size_t) totalChans);

            while (audioChannelsToUse.size() < totalChans)
                audioChannelsToUse.add (0);
        }

        template <typename FloatType>
        void perform (AudioBuffer<FloatType>& sharedBufferChans, const OwnedArray<MidiBuffer>& sharedMidiBuffers, const int numSamples)
        {
            auto& channels = audioChannels.get<FloatType>();

            for (int i = totalChans; --i >= 0;)
                channels[i] = sharedBufferChans.getWritePointer (audioChannelsToUse.getUnchecked (i), 0);

            AudioBuffer<FloatType> buffer (channels, totalChans, numSamples);

            if (processor->isSuspended())
            {
                buffer.clear();
            }
            else
            {
                ScopedLock lock (processor->getCallbackLock());

                callProcess (buffer, *sharedMidiBuffers.getUnchecked (midiBufferToUse));
            }
        }

        void callProcess (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
        {
            processor->processBlock (buffer, midiMessages);
        }

        void callProcess (AudioBuffer<double>& buffer, MidiBuffer& midiMessages)
        {
            if (processor->isUsingDoublePrecision())
            {
                processor->processBlock (buffer, midiMessages);
            }
            else
            {
                // if the processor is in single precision mode but the graph in double
                // precision then we need to convert between buffer formats. Note, that
                // this will only happen if the processor does not support double
                // precision processing.
                tempBuffer.makeCopyOf (buffer, true);
                processor->processBlock (tempBuffer, midiMessages);
                buffer.makeCopyOf (tempBuffer, true);
            }
        }

        const AudioProcessorGraph::Node::Ptr node;
        AudioProcessor* const processor;

    private:
        Array<int> audioChannelsToUse;
        FloatAndDoubleComposition<HeapBlock<FloatPlaceholder*>> audioChannels;
        AudioBuffer<float> tempBuffer;
        const int totalChans;
        const int midiBufferToUse;

        JUCE_DECLARE_NON_COPYABLE (ProcessBufferOp)
    };
};

//==============================================================================
template <typename RenderSequence>
struct SequenceBuilder
{
    SequenceBuilder (AudioProcessorGraph& g, RenderSequence& s)
        : graph (g), sequence (s)
    {
        createOrderedNodeList();

        nodeIds.add ((uint32) zeroNodeID); // first buffer is read-only zeros
        channels.add (0);

        midiNodeIds.add ((uint32) zeroNodeID);

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
    AudioProcessorGraph& graph;
    RenderSequence& sequence;

    Array<AudioProcessorGraph::Node*> orderedNodes;
    Array<int> channels;
    Array<uint32> nodeIds, midiNodeIds;

    enum { freeNodeID = 0xffffffff, zeroNodeID = 0xfffffffe, anonymousNodeID = 0xfffffffd };

    static bool isNodeBusy (uint32 nodeID) noexcept     { return nodeID != freeNodeID && nodeID != zeroNodeID; }

    Array<uint32> nodeDelayIDs;
    Array<int> nodeDelays;
    int totalLatency = 0;

    int getNodeDelay (uint32 nodeID) const        { return nodeDelays [nodeDelayIDs.indexOf (nodeID)]; }

    void setNodeDelay (uint32 nodeID, const int latency)
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

    int getInputLatencyForNode (uint32 nodeID) const
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

        bool isAnInputTo (uint32 sourceID, uint32 destID) const noexcept
        {
            return isAnInputToRecursive (sourceID, destID, entries.size());
        }

    private:
        //==============================================================================
        std::unordered_map<uint32, SortedSet<uint32>> entries;

        bool isAnInputToRecursive (uint32 sourceID, uint32 destID, size_t recursionDepth) const noexcept
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

        JUCE_DECLARE_NON_COPYABLE (NodeInputLookupTable)
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
            Array<uint32> sourceNodes;
            Array<int> sourceOutputChans;

            for (int i = graph.getNumConnections(); --i >= 0;)
            {
                auto* c = graph.getConnection (i);

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
                const uint32 srcNode = sourceNodes.getUnchecked(0);
                const int srcChan = sourceOutputChans.getUnchecked(0);

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

                    markBufferAsContaining (bufIndex, static_cast<uint32> (anonymousNodeID), 0);

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
        Array<uint32> midiSourceNodes;

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

        sequence.addProcessBufferOp (&node, audioChannelsToUse, totalChans, midiBufferToUse);
    }

    //==============================================================================
    int getFreeBuffer (const bool forMidi)
    {
        if (forMidi)
        {
            for (int i = 1; i < midiNodeIds.size(); ++i)
                if (midiNodeIds.getUnchecked(i) == freeNodeID)
                    return i;

            midiNodeIds.add ((uint32) freeNodeID);
            return midiNodeIds.size() - 1;
        }
        else
        {
            for (int i = 1; i < nodeIds.size(); ++i)
                if (nodeIds.getUnchecked(i) == freeNodeID)
                    return i;

            nodeIds.add ((uint32) freeNodeID);
            channels.add (0);
            return nodeIds.size() - 1;
        }
    }

    int getReadOnlyEmptyBuffer() const noexcept
    {
        return 0;
    }

    int getBufferContaining (const uint32 nodeId, const int outputChannel) const noexcept
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
                nodeIds.set (i, (uint32) freeNodeID);
            }
        }

        for (int i = 0; i < midiNodeIds.size(); ++i)
        {
            if (isNodeBusy (midiNodeIds.getUnchecked(i))
                 && ! isBufferNeededLater (stepIndex, -1,
                                           midiNodeIds.getUnchecked(i),
                                           AudioProcessorGraph::midiChannelIndex))
            {
                midiNodeIds.set (i, (uint32) freeNodeID);
            }
        }
    }

    bool isBufferNeededLater (int stepIndexToSearchFrom,
                              int inputChannelOfIndexToIgnore,
                              const uint32 nodeId,
                              const int outputChanIndex) const
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

    void markBufferAsContaining (int bufferNum, uint32 nodeId, int outputIndex)
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

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SequenceBuilder)
};

//==============================================================================
AudioProcessorGraph::Connection::Connection (const uint32 sourceID, const int sourceChannel,
                                             const uint32 destID, const int destChannel) noexcept
    : sourceNodeId (sourceID), sourceChannelIndex (sourceChannel),
      destNodeId (destID), destChannelIndex (destChannel)
{
}

//==============================================================================
AudioProcessorGraph::Node::Node (const uint32 nodeID, AudioProcessor* const p) noexcept
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
struct AudioProcessorGraph::AudioProcessorGraphBufferHelpers
{
    AudioProcessorGraphBufferHelpers()
    {
        currentAudioInputBuffer.floatVersion  = nullptr;
        currentAudioInputBuffer.doubleVersion = nullptr;
    }

    void setRenderingBufferSize (int newNumChannels, int newNumSamples)
    {
        renderingBuffers.floatVersion. setSize (newNumChannels, newNumSamples);
        renderingBuffers.doubleVersion.setSize (newNumChannels, newNumSamples);

        renderingBuffers.floatVersion. clear();
        renderingBuffers.doubleVersion.clear();
    }

    void release()
    {
        renderingBuffers.floatVersion. setSize (1, 1);
        renderingBuffers.doubleVersion.setSize (1, 1);

        currentAudioInputBuffer.floatVersion  = nullptr;
        currentAudioInputBuffer.doubleVersion = nullptr;

        currentAudioOutputBuffer.floatVersion. setSize (1, 1);
        currentAudioOutputBuffer.doubleVersion.setSize (1, 1);
    }

    void prepareInOutBuffers(int newNumChannels, int newNumSamples)
    {
        currentAudioInputBuffer.floatVersion  = nullptr;
        currentAudioInputBuffer.doubleVersion = nullptr;

        currentAudioOutputBuffer.floatVersion. setSize (newNumChannels, newNumSamples);
        currentAudioOutputBuffer.doubleVersion.setSize (newNumChannels, newNumSamples);
    }

    FloatAndDoubleComposition<AudioBuffer<FloatPlaceholder>> renderingBuffers;
    FloatAndDoubleComposition<AudioBuffer<FloatPlaceholder>*> currentAudioInputBuffer;
    FloatAndDoubleComposition<AudioBuffer<FloatPlaceholder>> currentAudioOutputBuffer;
};

//==============================================================================
AudioProcessorGraph::AudioProcessorGraph()
    : lastNodeId (0), audioBuffers (new AudioProcessorGraphBufferHelpers),
      currentMidiInputBuffer (nullptr), isPrepared (false)
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

AudioProcessorGraph::Node* AudioProcessorGraph::getNodeForId (const uint32 nodeId) const
{
    for (auto* n : nodes)
        if (n->nodeId == nodeId)
            return n;

    return nullptr;
}

AudioProcessorGraph::Node* AudioProcessorGraph::addNode (AudioProcessor* const newProcessor, uint32 nodeId)
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

bool AudioProcessorGraph::removeNode (const uint32 nodeId)
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
const AudioProcessorGraph::Connection* AudioProcessorGraph::getConnectionBetween (const uint32 sourceNodeId,
                                                                                  const int sourceChannelIndex,
                                                                                  const uint32 destNodeId,
                                                                                  const int destChannelIndex) const
{
    const Connection c (sourceNodeId, sourceChannelIndex, destNodeId, destChannelIndex);
    ConnectionSorter sorter;
    return connections [connections.indexOfSorted (sorter, &c)];
}

bool AudioProcessorGraph::isConnected (const uint32 possibleSourceNodeId,
                                       const uint32 possibleDestNodeId) const
{
    for (auto* c : connections)
        if (c->sourceNodeId == possibleSourceNodeId && c->destNodeId == possibleDestNodeId)
            return true;

    return false;
}

bool AudioProcessorGraph::canConnect (const uint32 sourceNodeId,
                                      const int sourceChannelIndex,
                                      const uint32 destNodeId,
                                      const int destChannelIndex) const
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

bool AudioProcessorGraph::addConnection (const uint32 sourceNodeId,
                                         const int sourceChannelIndex,
                                         const uint32 destNodeId,
                                         const int destChannelIndex)
{
    if (! canConnect (sourceNodeId, sourceChannelIndex, destNodeId, destChannelIndex))
        return false;

    ConnectionSorter sorter;
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

bool AudioProcessorGraph::removeConnection (const uint32 sourceNodeId, const int sourceChannelIndex,
                                            const uint32 destNodeId, const int destChannelIndex)
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

bool AudioProcessorGraph::disconnectNode (const uint32 nodeId)
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
    ScopedPointer<RenderSequence> sequence;

    {
        const ScopedLock sl (getCallbackLock());
        renderSequence.swapWith (sequence);
    }
}

bool AudioProcessorGraph::isAnInputTo (const uint32 possibleInputId,
                                       const uint32 possibleDestinationId,
                                       const int recursionCheck) const
{
    if (recursionCheck > 0)
    {
        for (auto* c : connections)
            if (c->destNodeId == possibleDestinationId
                 && (c->sourceNodeId == possibleInputId || isAnInputTo (possibleInputId, c->sourceNodeId, recursionCheck - 1)))
                return true;
    }

    return false;
}

void AudioProcessorGraph::buildRenderingSequence()
{
    ScopedPointer<RenderSequence> newRenderSequence;

    {
        MessageManagerLock mml;

        for (auto* node : nodes)
            node->prepare (getSampleRate(), getBlockSize(), this, getProcessingPrecision());

        newRenderSequence = new RenderSequence (*this);

        SequenceBuilder<RenderSequence> builder (*this, *newRenderSequence);
    }

    {
        // swap over to the new rendering sequence..
        const ScopedLock sl (getCallbackLock());

        audioBuffers->setRenderingBufferSize (newRenderSequence->numBuffersNeeded, getBlockSize());

        for (auto* b : midiBuffers)
            b->clear();

        while (midiBuffers.size() < newRenderSequence->numMidiBuffersNeeded)
            midiBuffers.add (new MidiBuffer());

        renderSequence.swapWith (newRenderSequence);
    }
}

void AudioProcessorGraph::handleAsyncUpdate()
{
    buildRenderingSequence();
}

//==============================================================================
void AudioProcessorGraph::prepareToPlay (double /*sampleRate*/, int estimatedSamplesPerBlock)
{
    audioBuffers->prepareInOutBuffers (jmax (1, getTotalNumOutputChannels()), estimatedSamplesPerBlock);

    currentMidiInputBuffer = nullptr;
    currentMidiOutputBuffer.clear();

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

    audioBuffers->release();
    midiBuffers.clear();

    currentMidiInputBuffer = nullptr;
    currentMidiOutputBuffer.clear();
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

template <typename FloatType>
void AudioProcessorGraph::processAudio (AudioBuffer<FloatType>& buffer, MidiBuffer& midiMessages)
{
    auto& renderingBuffers          = audioBuffers->renderingBuffers.get<FloatType>();
    auto& currentAudioInputBuffer   = audioBuffers->currentAudioInputBuffer.get<FloatType>();
    auto& currentAudioOutputBuffer  = audioBuffers->currentAudioOutputBuffer.get<FloatType>();

    const int numSamples = buffer.getNumSamples();
    jassert (numSamples <= getBlockSize());

    currentAudioInputBuffer = &buffer;
    currentAudioOutputBuffer.setSize (jmax (1, buffer.getNumChannels()), numSamples);
    currentAudioOutputBuffer.clear();
    currentMidiInputBuffer = &midiMessages;
    currentMidiOutputBuffer.clear();

    if (renderSequence != nullptr)
        renderSequence->perform (renderingBuffers, midiBuffers, numSamples);

    for (int i = 0; i < buffer.getNumChannels(); ++i)
        buffer.copyFrom (i, 0, currentAudioOutputBuffer, i, 0, numSamples);

    midiMessages.clear();
    midiMessages.addEvents (currentMidiOutputBuffer, 0, buffer.getNumSamples(), 0);
}

template <typename FloatType>
void AudioProcessorGraph::sliceAndProcess (AudioBuffer<FloatType>& buffer, MidiBuffer& midiMessages)
{
    auto n = buffer.getNumSamples();
    auto ch = buffer.getNumChannels();
    auto max = 0;

    for (auto pos = 0; pos < n; pos += max)
    {
        max = jmin (n - pos, getBlockSize());

        AudioBuffer<FloatType> audioSlice (buffer.getArrayOfWritePointers(), ch, pos, max);
        MidiBuffer midiSlice;

        midiSlice.addEvents (midiMessages, pos, max, 0);
        processAudio (audioSlice, midiSlice);
    }
}

double AudioProcessorGraph::getTailLengthSeconds() const            { return 0; }
bool AudioProcessorGraph::acceptsMidi() const                       { return true; }
bool AudioProcessorGraph::producesMidi() const                      { return true; }
void AudioProcessorGraph::getStateInformation (juce::MemoryBlock&)  {}
void AudioProcessorGraph::setStateInformation (const void*, int)    {}

void AudioProcessorGraph::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    sliceAndProcess (buffer, midiMessages);
}

void AudioProcessorGraph::processBlock (AudioBuffer<double>& buffer, MidiBuffer& midiMessages)
{
    sliceAndProcess (buffer, midiMessages);
}

// explicit template instantiation
template void AudioProcessorGraph::processAudio<float> ( AudioBuffer<float>& buffer,
                                                         MidiBuffer& midiMessages);
template void AudioProcessorGraph::processAudio<double> (AudioBuffer<double>& buffer,
                                                         MidiBuffer& midiMessages);

//==============================================================================
AudioProcessorGraph::AudioGraphIOProcessor::AudioGraphIOProcessor (const IODeviceType deviceType)
    : type (deviceType), graph (nullptr)
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

template <typename FloatType>
void AudioProcessorGraph::AudioGraphIOProcessor::processAudio (AudioBuffer<FloatType>& buffer,
                                                               MidiBuffer& midiMessages)
{
    auto& currentAudioInputBuffer  = graph->audioBuffers->currentAudioInputBuffer.get<FloatType>();
    auto& currentAudioOutputBuffer = graph->audioBuffers->currentAudioOutputBuffer.get<FloatType>();

    jassert (graph != nullptr);

    switch (type)
    {
        case audioOutputNode:
        {
            for (int i = jmin (currentAudioOutputBuffer.getNumChannels(),
                               buffer.getNumChannels()); --i >= 0;)
            {
                currentAudioOutputBuffer.addFrom (i, 0, buffer, i, 0, buffer.getNumSamples());
            }

            break;
        }

        case audioInputNode:
        {
            for (int i = jmin (currentAudioInputBuffer->getNumChannels(),
                               buffer.getNumChannels()); --i >= 0;)
            {
                buffer.copyFrom (i, 0, *currentAudioInputBuffer, i, 0, buffer.getNumSamples());
            }

            break;
        }

        case midiOutputNode:
            graph->currentMidiOutputBuffer.addEvents (midiMessages, 0, buffer.getNumSamples(), 0);
            break;

        case midiInputNode:
            midiMessages.addEvents (*graph->currentMidiInputBuffer, 0, buffer.getNumSamples(), 0);
            break;

        default:
            break;
    }
}

void AudioProcessorGraph::AudioGraphIOProcessor::processBlock (AudioBuffer<float>& buffer,
                                                               MidiBuffer& midiMessages)
{
    processAudio (buffer, midiMessages);
}

void AudioProcessorGraph::AudioGraphIOProcessor::processBlock (AudioBuffer<double>& buffer,
                                                               MidiBuffer& midiMessages)
{
    processAudio (buffer, midiMessages);
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
