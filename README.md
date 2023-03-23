# seq_geom_xform

This crate consists of both a library and program to aid in processing
sequencing libraries with different "fragment geometries". This crate is very
much under active development, and so suggestions and feedback are welcome.
While we do aim to make this library and tool as general as possible,
development is driven by our primary use case in supporting the most common
geometries present in single-cell sequencing data.

The goal of this crate, is to consume a sequencing library and 
a [sequence fragment geometry description specification](https://hackmd.io/@PI7Og0l1ReeBZu_pjQGUQQ/rJMgmvr13)
and to then parse the library in accordance with the description.  Specifically, 
this tool is most useful when one has a "complex" geometry (i.e. a geometry where the 
position or length of some sequence segment — a UMI or cellular barcode — is not fixed) 
and they need to transform the library into one that encodes equivalent information 
int a "simple" geometry (i.e. a geometry where all sequence segments are at fixed 
and known positions and of a precisely known length).

# Basic usage

```
Transform/normalize complex single-cell fragment geometries into simple geometries.

Usage: seq_xformer [OPTIONS] --geom <GEOM> --out1 <OUT1> --out2 <OUT2>

Options:
  -g, --geom <GEOM>    Expected input read geometry specification
  -1, --read1 <READ1>  read 1 files, comma delimited
  -2, --read2 <READ2>  read 2 files, comma delimited
  -o, --out1 <OUT1>    where output r1 should be written (currently uncompressed)
  -w, --out2 <OUT2>    where output r2 should be written (currently uncompressed)
  -h, --help           Print help
  -V, --version        Print version
```

The `seq_xformer` program takes as input a [sequence fragment geometry
description specification](https://hackmd.io/@PI7Og0l1ReeBZu_pjQGUQQ/rJMgmvr13)
and a pair of input libraries (i.e. one or more paired-end files).  It will
then write the transformed sequences to the specified output files `--out1` and
`--out2`.  These could be regular files on disk, or, if you wish, they could be
[`fifos`](https://www.ibm.com/docs/en/aix/7.1?topic=m-mkfifo-command) that you
have set up for some receiving program to read from. The `seq_xformer` tool works 
in a streaming fashion, and so read pairs will be read from the input, transformed
and directly written to the output.


## Normalization

The normalization of complex geometries in the context of `seq_xformer` consists of 
turning variable-length segments into fixed-length segments, determining the position 
and content of variable position segments (most often determined by an anchor sequence),
and outputting a transformed sequence where every sequence segment is at a fixed and 
known position and has a single, fixed length.  Additionally, non-functional sequence 
components (e.g. anchor sequences or other padding) is removed.

## Transformation of variable-length segments

If an input library contains a variable-length segment (e.g. a segment that has
a minimum and maximum possible length that differ), then `seq_xformer` has a
specific strategy for turning these into fixed-length segements.  Specifically,
this is done by padding variable length segments so that no padded segments of
different lengths will collide.  For example, suppose that we have a segment
that constitutes part of a cellular barcode, and that this segment is of some
length between 8 and 10.  That is, when we see this segment, it will always have
length at least 8, and it will never have length more than 10.

In this case, `seq_xformer` will transform this variable length segment in the 
input into a fixed length segment of length 11 (the maximum length + 1) in the 
output.  This is done with the following padding strategy.  If an observed
segment in the input is of the maximum length (here 10), an `A` is appended
to it before it is written in the output.  If an observed segment is of length 
9, then `AC` is appended to it.  If an observed segment is of length 8, then 
`AAG` is appended to it.  Here, you can see that, since these segments are 
all padded with nucleotide strings of various length, they all end up having 
a fixed length (in this case 11) in the output.  Further, because observed
 input segments of every distinct length have a padding sequence that differs
 in the last character, then segments with a different initial lengths, 
 by construction, cannot collide.  Currently, `seq_xformer` supports segments
 whose length varies by up to 4 bases.  The general strategy is more scalable 
 (i.e. if the output length was the maximum input length + 2, then the variable 
 length window could be doubled, etc.).  However, this length restriction is 
 only enforced "per-piece".  So, for example, if a cellular barcode was split 
 across 2 separate segments, then each could have a length that varies by 
 up to 4 nucleotides.
