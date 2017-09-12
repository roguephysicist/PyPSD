# Usage: python lines.py
# then run gnuplot on the individual gnuplot files to generate the plots in PDF.

import numpy as np

params = [('PAN500KD_1', 17, 2), # done 000
          ('PAN500KD_2',  5, 2), # done 0.86
          ('PES200KD_1',  7, 3), # done 0.77
          ('PES200KD_2',  9, 3), # done 1.51
          ('PES300KD_1',  2, 3), # done 0.37
          ('PES300KD_2',  6, 2), # done 1.22
           ('PVDF04u_1',  2, 4), # done 0.42
           ('PVDF04u_2',  3, 3), # done 0.66
           ('PVDF15u_1',  5, 3), # done 0.57
           ('PVDF15u_2',  2, 2)] # done 0.52

for case in params:
    file = case[0]
    data = np.loadtxt('data/' + file + '.txt')
    elem = len(data)-1
    
    start_data = np.take(data, list(range(case[1])), 0)
    start_line = np.polyfit(start_data[:,0],start_data[:,2],1)
    
    end_data = np.take(data, list(range(elem, elem-case[2], -1)), 0)
    end_line = np.polyfit(end_data[:,0],end_data[:,2],1)
    
    intersect_x = (end_line[1] - start_line[1])/\
                  (start_line[0] - end_line[0])
    intersect_y = (start_line[0]*end_line[1] - end_line[0]*start_line[1])/\
                  (start_line[0] - end_line[0])
    intersection = [intersect_x, intersect_y]
    
    with open('gnuplot_' + str(file) + '.gp', 'w') as outfile:
        print('set terminal pdfcairo  transparent enhanced fontscale 0.5 size 5.00in, 3.00in', file=outfile)
        print('set output "{0}.pdf"'.format(file), file=outfile)
        # print('set terminal pngcairo  background "#ffffff" enhanced fontscale 1.0 size 1200, 720', file=outfile)
        # print('set output "{0}.png"'.format(file), file=outfile)
        # print('set locale "en_US.UTF-8"', file=outfile)
        print('GNUTERM = "wxt"', file=outfile)
        print('set grid', file=outfile)
        print('set xrange [0:*]', file=outfile)                        
        print('set yrange [0:*]', file=outfile)                        
        print('set xlabel "TMP (bar)"', file=outfile)                        
        print('set ylabel "J (L/h/m^{2})"', file=outfile)
        print('f(x) = (x < ({0:.4f} + {0:.4f}*0.1) ) ? {1:.4f}*x + {2:.4f} : 1/0'.format(intersection[0], start_line[0], start_line[1]), file=outfile)
        print('g(x) = (x > ({0:.4f} - {0:.4f}*0.1) ) ? {1:.4f}*x + {2:.4f} : 1/0'.format(intersection[0], end_line[0], end_line[1]), file=outfile)
        print('set label at {0:.4f},{1:.4f} "" point lw 2 pt 6 ps 1 front'.format(intersection[0], intersection[1]), file=outfile)
        print('set label at {0:.4f},0 "{0:.2f}" offset 1,2 front'.format(intersection[0], intersection[1]), file=outfile)
        print('set arrow from {0:.4f},{1:.4f} to {0:.4f},0 nohead lw 2 dt 2 front'.format(intersection[0], intersection[1]), file=outfile)
        print('p "data/{0}.txt" u 1:3:4 t "" w errorbars lw 1.5 pt 7 ps 0.5 lc rgb "#268bd2",\\'.format(file), file=outfile)
        # print('"" u 1:3 t "" lw 1.5 lc rgb "#268bd2" w l,\\'.format(file), file=outfile)
        print('f(x) t "" lw 1.5 lc rgb "#dc322f" w l,\\', file=outfile)
        print('g(x) t "" lw 1.5 lc rgb "#dc322f" w l', file=outfile)
