import os

def check_if_file_exists(name):
    ret_val = False
    if os.path.isfile(name):
        ret_val = True
    return ret_val
    
def string_to_chars(word):
    return [char for char in word]

def madx_plot_file(seqname, filename):
    l1 = 'setplot, font=4, xsize=34, ysize=25;\n'
    l2 = "plot, table=twiss, haxis=s, vaxis1=betx,bety,dx,VAXIS2=MUX,MUY,INTERPOLATE,vmin=0,vmax=30.,title=sequence:"+seqname+", style=100, FILE="+filename+';\n' 
    plot_command = l1 + l2
    return plot_command
    
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def madx_plot_file(seqname, filename):
    l1 = 'setplot, font=4, xsize=34, ysize=25;\n'
    l2 = "plot, table=twiss, haxis=s, vaxis1=betx,bety,dx,vaxis2=mux,muy,interpolate,vmin=0,vmax=30.,title=sequence:"+seqname+", style=100, FILE="+filename+';\n' 
    plot_command = l1 + l2
    return plot_command

def convert_mad8_to_madx(inpath, mad8file, outpath=None, madxfile=None):

    if madxfile is None:  madxfile = mad8file + '.madx'
    if outpath is None:   outpath = inpath
    if outpath is inpath and madxfile is mad8file:
        print('ERROR: MAD-8 to MAD-X Conversion: input and output file have the same name: Exiting')
        exit(0)
        
    if not os.path.isdir(outpath): os.makedirs(outpath)        
        
    mad8file = inpath + mad8file
    madxfile = outpath + madxfile
    
    if check_if_file_exists(madxfile):
        print(madxfile, ' already exists, deleting')
        os.remove(madxfile)
    f8 = open(mad8file, 'r')
    lines = []
    
    infile = open(mad8file)
    linelist = infile.readlines()
    file_len = len(linelist)
    
    skip_lines = []
    last_sequence = ''
    twiss_files = []
    
    # Using readlines we can manually iterate over the lines and thus skip lines at need
    for i in range(0,file_len,1):
        if i in skip_lines:
            pass
        else:
            line = linelist[i]
            
            if line.startswith('!'):
                #line += '\n'
                lines.append(line)
                pass
            
            else:
                if not line.isspace():
                    # Handle calling of other files - case sensitive - don't use .lower() if it is a file call as filenames are case sensitive
                    if 'call,' in line:
                        line = line.strip()                    
                        line = line.replace('call,', 'call,file=')
                    elif 'CALL,' in line:
                        line = line.strip()                    
                        line = line.replace('CALL,', 'CALL,FILE=')
                    else:
                        line = line.lower().strip()            

                    # SPECIAL CASE - MULTIPOLE definition
                    # Note that MAD8 uses KnL with tilts as Tn
                    # MADX uses KnL and KnS for normal and skew components
                    # Here we assume no tilts in MAD8 (thus no skew components in MAD-X)
                    if 'multipole' in line:
                        N_poles = []
                        N_pole_val = []
                        #S_poles = []
                        #S_pole_val = []
                        line_multi = ''

                        # first check if we have a multiple line definition
                        if line.endswith('&'):
                            # Reduce all lines in the multipole definition to one for processing

                            #print('multipole')
                            multi_def = True
                            while multi_def:
                                if '&' not in line : multi_def = False
                                line = line.strip('\n')
                                line = line.strip('\t')
                                line = line.strip(' ')
                                line = line.strip('&')
                                line_multi = line_multi + line
                                skip_lines.append(i)
                                i+=1
                                line = linelist[i]
                                line = line.lower().strip()
                            line_multi = line_multi.strip('&')
                            line_multi = line_multi.strip('\n')
                            line_multi = line_multi.strip('\t')
                            line_multi = line_multi.strip(' ')
                            line_multi = line_multi.replace(' ', '')
                            line_multi = line_multi.replace('\t', '')
                            #print(line_multi)           

                        multiline = line_multi.split(',')
                        multi_name = multiline[0].replace(':', ' : ')
                        #print(multi_name)                
                        madx_multi_line = multi_name

                        for ml in multiline:
                            if '=' in ml:
                                multi_def = ml.split('=')
                                #print(multi_def)

                                if 'l' in multi_def[0]:
                                    for s in string_to_chars(multi_def[0]): 
                                        if s.isdigit():
                                            N_poles.append(int(s))
                                    N_pole_val.append(str(multi_def[1]))
                                #elif 's' in multi_def[0]:
                                    #for s in string_to_chars(multi_def[0]): 
                                        #if s.isdigit():
                                            #S_poles.append(int(s))
                                    #S_pole_val.append(str(multi_def[1]))
                                else:
                                    print('ERROR: MAD-8 to MAD-X Conversion: Multipole not as expected')
                                    print(line)
                                    exit(0)

                        #print(N_poles)
                        #print(N_pole_val)

                        madx_multi_line = madx_multi_line + ', knl={'

                        N_pole_it = 0
                        for i in range(0,11,1):
                            if i in N_poles:
                                madx_multi_line = madx_multi_line + N_pole_val[N_pole_it] + ', '
                                N_pole_it+=1
                            else:
                                madx_multi_line = madx_multi_line + '0, '

                        madx_multi_line = madx_multi_line.strip()[:-1]
                        madx_multi_line += '}'
                        lines.append(madx_multi_line)         

                    
                    # SPECIAL CASE - VARY
                    if 'vary' in line:
                        line = line.replace('vary,', 'vary, name=')
                    
                    # SPECIAL CASE - CONSTRAINT - may need manual inspection if range/seq not obvious
                    if 'constraint' in line:
                        if '#e' in line.split(',')[1]:
                            line = line.replace('#e', 'range=#e')
                        else:
                            line = line.replace('constraint,', 'constraint, sequence=')                            
                                               
                    # Handle line endings
                    if line.endswith('&'):
                        line = line[:-1]
                    else:
                        if '!' in line: # line with a comment in - assume only one !
                            commented_line = line.split('!')
                            line = commented_line[0].strip() + '; !' + commented_line[1] + '\n'
                        else:
                            line += ';\n'
                    
                    # SPECIAL CASE - USE - store last sequence for twiss file        
                    if 'use,' in line:
                        last_sequence = line.split(',')[1].replace(';','').replace('\n','')
                        print ('Last sequence is ', last_sequence)                        

                    # SPECIAL CASE - TWISS file output
                    if 'twiss,save' in line:                        
                        twiss_counter = 0
                        madx_string = 'select, flag=twiss, clear;\nselect, flag=twiss, column=name, keyword, s, l, mux, muy, betx, alfx, bety, alfy, x, px, y, py, t, pt, dx, dpx, dy, dpy, energy, angle, k0l, k0sl, k1l, k1sl, k2l, k2sl, k3l, k3sl, k4l, k4sl, k5l, k5sl, vkick, hkick, tilt, slot_id, volt, lag, freq, harmon, beta11, beta12, alfa11, alfa12, disp1, disp2, disp3, disp4;\nset, format="12.12f";\n'
                        twiss_file = 'twiss_' + last_sequence + '_' + str(twiss_counter) +'.tfs'
                        if twiss_file in twiss_files:
                            while twiss_file in twiss_files:
                                twiss_counter+=1
                                twiss_file = 'twiss_' + last_sequence + '_' + str(twiss_counter) +'.tfs'
                        else: twiss_files.append(twiss_file)
                                
                        madx_end_string = 'twiss, sequence=' + last_sequence + ', file=' + twiss_file + ', save;\n'
                                                    
                        line = madx_string + madx_end_string                        

                    # SPECIAL CASE - PLOT
                    if 'plot' in line: 
                        if line.endswith(','): 
                            line = line + linelist[i+1]
                            skip_lines.append(i+1)
                        line = madx_plot_file(last_sequence,'plot')    
                        
                    # Change [] access to ->
                    if '[' and ']' in line:
                        line = line.replace('[', '->')
                        line = line.replace(']', '')

                    line = line.replace('quad', 'quadrupole')
                    line = line.replace('line', 'line=')
                    line = line.replace('spline=', 'interpolate')
                    line = line.replace('use,', 'use, sequence=')
                    line = line.replace('filename', 'file')
                    line = line.replace('mad', outpath + 'mad')
                    
                    if line.startswith("print"): line = '!' + line
                    if line.startswith("assign,echo"): line = '!' + line
                    if 'optics' in line: line = '!' + line
                    if 'survey' and 'tape=' in line: line = '!' + line
                        
                    # Element checks - add ':=' when using variables as values
                    if 'sbend' in line:
                        if ' l' in line:
                            if not is_number(line.split(' l')[1].split('=')[1].split(',')[0]):
                                line = line.replace(' l=', ' l:=')
                                
                        if 'k1' in line:
                            if not is_number(line.split('k1')[1].split('=')[1].split(',')[0]):
                                line = line.replace('k1=', 'k1:=')
                        
                        if 'angle' in line:
                            if not is_number(line.split('angle')[1].split('=')[1].split(',')[0]):
                                line = line.replace('angle=', 'angle:=')
                        print(':= inserted for variable:', line)
                        
                        
                    if 'quadrupole' in line:
                        if not is_number(line.split('k1')[1].split('=')[1].split(';')[0]):
                            line = line.replace('k1=', 'k1:=')
                        print(':= inserted for variable:', line)
                        
                    if 'sextupole' in line:
                        if not is_number(line.split('k2')[1].split('=')[1].split(';')[0]):
                            line = line.replace('k2=', 'k2:=')
                        print(':= inserted for variable:', line)
                        
                    if 'kicker' in line:
                        if not is_number(line.split('kick')[2].split('=')[1].split(';')[0]):
                            line = line.replace('kick=', 'kick:=')
                        print(':= inserted for variable:', line)                                           
                                                
                lines.append(line)
            
    with open(madxfile, 'x') as fx:
        for line in lines:
            fx.write(line)

# Example of how to use the function

madfile = 'ISIS_II_EHRCS.mad8'
madfile2 = 'ISIS_II_EHRCS.madx'
elfile = 'ISIS_II_EHRCS.elements'
strfile = 'ISIS_II_EHRCS.strength'
seqfile = 'ISIS_II_EHRCS.sequence'
beamfile = 'ISIS_II_EHRCS.beam'
inpath = '../00_Lattice_Files/Master_Lattice_Files/MAD8/'
outpath = '../00_Lattice_Files/Master_Lattice_Files/MADX/'

convert_mad8_to_madx(inpath, madfile, outpath, madfile2)
convert_mad8_to_madx(inpath, elfile, outpath, elfile)
convert_mad8_to_madx(inpath, strfile, outpath, strfile)
convert_mad8_to_madx(inpath, seqfile, outpath, seqfile)
convert_mad8_to_madx(inpath, beamfile, outpath, beamfile)
