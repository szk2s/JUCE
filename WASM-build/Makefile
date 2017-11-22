# Web Audio Modules (WAMs)
# jari kleimola 2017 (jari@webaudiomodules.org)
# wasm makefile for obxd

TARGET = ./obxd.js
SRC = ../../Source
WAM = ./wamsdk
JUCE = ../../JuceLibraryCode
INC = -I. -I$(SRC) -I$(JUCE)  -I$(JUCE)/modules

SOURCES = \
$(SRC)/PluginProcessor.cpp \
$(WAM)/processor.cpp \
$(JUCE)/modules/juce_audio_plugin_client/WAM/juce_WAM_Wrapper.cpp

LIBS = lib/juce-audioworklet.bc

Wno = -Wno-inconsistent-missing-override -Wno-deprecated-register -Wno-unused-private-field \
-Wno-reorder -Wno-implicit-conversion-floating-point-to-bool -Wno-unused-variable

JUCEFLAGS = -DJUCE_STANDALONE_APPLICATION=0 -DJUCE_GLOBAL_MODULE_SETTINGS_INCLUDED=1
CFLAGS  = -DAUDIOWORKLET -std=c++11 -Wall $(Wno) $(INC) $(JUCEFLAGS)
LDFLAGS = -O2
EMFLAGS = -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXPORT_NAME="'AudioWorkletGlobalScope.WAM'"

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(EMFLAGS) $(LIBS) -o $@ $(SOURCES)
