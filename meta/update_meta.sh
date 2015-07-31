#!/bin/sh
# enums
g++ -o enum_translator main_enum_translator.cpp
./enum_translator
rm enum_translator
# objects
g++ -o object_translator main_object_translator.cpp
./object_translator
rm object_translator

