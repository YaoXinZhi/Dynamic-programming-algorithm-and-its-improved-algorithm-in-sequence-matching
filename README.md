This folder contains the codes of the three algorithms, the code instructions README.md file, and project reports.

The correspondence between the three codes is

align_DP.py corresponds to the original dynamic programming algorithm.
align_band.py is an improved algorithm for dynamic programming of band-DP.
align_drop.py corresponds to X-drop's improved dynamic programming algorithm.

You can run the code directly through the command line:

> python align_DP.py parameter.txt input.txt output.txt
> python align_band.py parameter.txt input.txt output.txt
> python align_drop.py parameter.txt input.txt output.txt

The parameter file specifies the parameters of the three algorithms, and the penalty matrix.
The parameter file format is as follows:

1; the threshold X for X-drop
3; bandwidth B
0; score for initiating a gap
-1; score for each base insert/delete

; Below is the alphabet
a c g t

; Below is the similarity matrix for the alphabet
2 -1 -1 -1
-1 2 -1 -1
-1 -1 2 -1
-1 -1 -1 2

The input file contains the two sequences to be compared.
The input file format is as follows:

; This is the example in slide 45
>seq1
tgacaatccc

>seq2
tgagcatggt


The output file is the result save file.
The output file format is as follows:

score = 10
entries = 121

>seq1

tga-caat

>seq2

tgagc-at

Where score is the alignment score, entries are the effective calculation positions, and seq1 and seq2 correspond to the two sequences to be aligned, respectively.

If you have any questions, please contact my mailbox, xinzhi_bioinfo@163.com, or visit my Github, https://github.com/YaoXinZhi.

Best Wishes.


