# Tracés gnuplot avec styles de lignes perso

# Output dans un fichier pdf
set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 10,6 # inches

#set term postscript eps enhanced color linewidth 2 font "Arial, 12"
#    size 10,6 # inches

# Styles de ligne pour les axes
set style line 100 lw 2 lt rgb "#808080" # grey

# Style de ligne pour la grille
set style line 101 lt 0             # dashed
set style line 101 lt rgb "#808080" # grey

set grid back linestyle 101

set border 3 back linestyle 100 # Remove border on top and right.  These
             # borders are useless and make it harder
             # to see plotted lines near the border.
	     # Also, put it in grey; no need for so much emphasis on a border.

set xtics nomirror
set ytics nomirror



# Line styles: try to pick pleasing colors, rather
# than strictly primary colors or hard-to-see colors
# like gnuplot's default yellow.  Make the lines thick
# so they're easy to see in small plots in papers.
set style line 1 lt rgb "#DB0000" lw 2 # m Rouge
set style line 2 lt rgb "#00009E" lw 2 # f Bleu foncé
set style line 3 lt rgb "#00994D" lw 2 # x1 Vert
set style line 4 lt rgb "#F25900" lw 2 # h dashed Orange
set style line 5 lt rgb "#DB006E" lw 2 # j Rouge rose
set style line 6 lt rgb "#7ACAFF" lw 2 # d Bleu clair
set style line 7 lt rgb "#FF66FF" lw 2 # Rose
set style line 10 linetype 1 linecolor rgb "#00994D" lw 2 #2
set style line 11 linetype 1 linecolor rgb  "#DB0000"  lw 1
set style line 12 linetype 1 linecolor rgb "#5080cb" lw 3
set style line 13 linetype 1 linecolor rgb  "#FF66FF"  lw 3


#set title "sin"
set key bottom right

set xlabel "x(m)"
set ylabel "y(m)"

i = 0

set output 'SOL1D_proc0.0.pdf'

plot   'SOL1D_proc0.0.res'  u 1:5  w  lp  ls 10 t 'Surface eau' ,  'SOL1D_proc0.0.res'  u 1:4  w  lp  ls 11 t 'Fond' 

i = 1000

fil(n) = sprintf('SOL1D_proc0.%04.0f.res',n)

ou(n) = sprintf('SOL1D_proc0.%04.0f.pdf',n)

while (i <= 9999 ){

set output ou(i)

plot  fil(i) u 1:5  w  lp  ls 10 t 'Surface eau' ,  fil(i) u 1:4  w  lp  ls 11 t 'Fond'
i = i + 1000
}