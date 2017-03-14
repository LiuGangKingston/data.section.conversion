# data.section.conversion

This is a template code for converting some identified sections
of data in a text file and keep the rest part unchanged for
general purposes. In other words, an output file will be generated
with some data converted and all the rest simply copied from
the input file. Specifically, this code converts Cartesian
coordinates into fractionals coordinates of particles in a cell
with respect to the three period vectors of a periodic-structured
system, as an example. Although not necessarily, it is assumed that
each section to be converted contains all data needed for the
converting operation, which means in this example all the period
vectors are also available in each section. The identifier of each
section should be in the first line of the section, but no data to
be converted in that line. Sections of the same structure/format
and the same way for conversion can be identified the same, then
converted by the same routine(s) of the code, supplied by the user.
Those sections are regarded as in the same group. Different sections
of the same group may have different amount of data to be converted.
If the data section does not have a unique identifier in the original
input file, a line with such can be inserted into the input file,
then the user has a choice to copy such a line into the output file
or not.

To use this template code, in most cases, users should edit

     1) the KEY_WORD_GROUPS module;

     2) SUBROUTINE DIRECTION();

     3) SUBROUTINEs for specific conversion of each group of data sections.


The KEY_WORD_GROUPS module is used to specify total number
of groups of data sections to be converted (3 here):

    the key words (identifier) for each group;

    whether the key word line will be copied into the output file;

    the number of lines after the key word line should be copied into the output file directly, 
    which do not contain any data to be converted and/or needed.


Following the example routines at the end

    SUBROUTINE NICE_MOLECULE_PROCESSING()

    SUBROUTINE XYZ_FORMAT

    SUBROUTINE GEOMETRY()

 users should code their own routines of such to perform their specific and detailed conversion. For that purposes, many tools can be used from the FOR_READ_AND_WRITE module. Once they are coded, they should be registered in the routione

    SUBROUTINE DIRECTION()

 based on the groups they are dealing with.



Any comments please send to gang.liu@queensu.ca

Copyright (c) CAC (HPCVL), Queen's University, Mar. 2017

