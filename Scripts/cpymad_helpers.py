import os
import math
import subprocess
import numpy as np
import pandas as pnd

from cpymad.madx import Madx
from cpymad.madx import Sequence
from cpymad.madx import SequenceMap
import tfs

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

import helper_functions

########################################################################
# lowercase column names for a given dataframe
########################################################################
def pandas_dataframe_lowercase_columns(dataframe):
    dataframe.columns = map(str.lower, dataframe.columns)

########################################################################
# write text to cpymad logfile
########################################################################
def cpymad_write_to_logfile(cpymad_logfile, log_string):
    f = open(cpymad_logfile, 'a')
    f.write('\n')
    f.write(log_string)
    f.close()
    
########################################################################
# start cpymad run with output to logfile
########################################################################
def cpymad_start(cpymad_logfile = './cpymad_logfile.log'):
    f = open(cpymad_logfile, 'w')
    madx_instance = Madx(stdout=f)      
    madx_instance.options.echo=True
    madx_instance.options.warn=True
    
    log_string = '! cpymad_start called'
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    return madx_instance

########################################################################
# print cpymad logfile text
########################################################################
def cpymad_print_output(cpymad_logfile = './cpymad_logfile.log', from_line = None):
    f = open(cpymad_logfile, "r")
    j = 0
    
    if from_line is None:
        from_line = 0
    else:
        j = from_line
        
    file_lines = f.readlines()
    
    for i in enumerate(file_lines):
        if int(i[0]) >= int(from_line):
            line = file_lines[i[0]].replace('\n','')
            print(line)
            j += 1
  
    return j

########################################################################
# return cpymad logfile text as list of lines (character strings)
########################################################################    
def cpymad_get_output(cpymad_logfile = './cpymad_logfile.log', from_line = None):
    f = open(cpymad_logfile, "r")
    j = 0
    
    if from_line is None:
        from_line = 0
    else:
        j = from_line
        
    file_lines = f.readlines()
    final_file_lines = []
    
    for i in enumerate(file_lines):
        if int(i[0]) >= int(from_line):
            line = file_lines[i[0]].replace('\n','')
            final_file_lines.append(line)
            j += 1
  
    return final_file_lines, j
 
########################################################################
# print active sequence
########################################################################       
def cpymad_get_active_sequence(madx_instance): return SequenceMap(madx_instance)

########################################################################
# return active sequence (first in list) as string
########################################################################     
def cpymad_get_active_sequence_name(madx_instance): return str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]

########################################################################
# return sequence list, and dict of elements
########################################################################     
def cpymad_sequence_to_dict(madx_instance, cpymad_logfile, madx_sequence):
    elements = dict()
    sequence = []
    for el in madx_sequence.elements:
        sequence.append(el.name)
        if el.name not in elements:
            elements[el.name] = el      
    
    log_string = '! cpymad_sequence_to_dict called for sequence '+ str(madx_sequence.name)
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    return sequence,elements

########################################################################
# check if sequence exists
########################################################################    
def cpymad_check_and_use_sequence(madx_instance, cpymad_logfile, sequence_name):     
    madx_instance.use(sequence=sequence_name)
    if 'warning' and sequence_name in cpymad_get_output(cpymad_logfile)[0][-1]:
        print(cpymad_get_output(cpymad_logfile)[0][-1])
        print('cpymad_check_and_use_sequence::Sequence not valid in this instance of MAD-X')           
        log_string = '! cpymad_check_and_use_sequence called for sequence ' + sequence_name
        cpymad_write_to_logfile(cpymad_logfile, log_string)  
        return False
    else: 
        print('Sequence',sequence_name,'exists in this instance of MAD-X\nActive sequence:\n')
        print(cpymad_get_active_sequence(madx_instance))       
        log_string = '! cpymad_check_and_use_sequence called for sequence ' + sequence_name
        cpymad_write_to_logfile(cpymad_logfile, log_string)  
        return True
    

########################################################################
# perform madx twiss with use sequence, return twiss dataframe
########################################################################    
def cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name, file_out=None):  
    if cpymad_check_and_use_sequence(madx_instance, cpymad_logfile, sequence_name):
        
        log_string = '! cpymad_madx_twiss called for sequence ' + sequence_name
        cpymad_write_to_logfile(cpymad_logfile, log_string)
        
        if file_out is None: file_out = sequence_name +'_madx_twiss.tfs'
    
        madx_instance.input('set, format="12.12f"')
        madx_instance.input('select, flag=twiss, column=keyword, name, s, l, betx, alfx, mux, bety, alfy, muy, x, px, y, py, t, pt, dx, dpx, dy, dpy, wx, phix, dmux, wy, phiy, dmuy, ddx, ddpx, ddy, ddpy, r11, r12, r21, r22, energy, angle, k0l, k0sl, k1l, k1sl, k2l, k2sl, k3l, k3sl, k4l, k4sl, k5l, k5sl, k6l, k6sl, k7l, k7sl, k8l, k8sl, k9l, k9sl, k10l, k10sl, ksi, hkick, vkick, tilt, e1, e2, h1, h2, hgap, fint, fintx, volt, lag, freq, harmon, slot_id, assembly_id, mech_sep, kmax, kmin, calib, polarity, alfa, beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33, alfa11, alfa12, alfa13, alfa21, alfa22, disp1, disp2, disp3, disp4')
        madx_instance.twiss(sequence=sequence_name, file=file_out)
        
        return madx_instance.table.twiss.dframe()
        
########################################################################
# perform madx twiss without use sequence, return twiss dataframe
########################################################################  
def cpymad_madx_twiss_nocheck(madx_instance, cpymad_logfile, sequence_name, file_out=None):  
    log_string = '! cpymad_madx_twiss_nocheck called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)

    if file_out is None: file_out = sequence_name +'_madx_twiss.tfs'

    madx_instance.input('set, format="12.12f"')
    madx_instance.input('select, flag=twiss, column=keyword, name, s, l, betx, alfx, mux, bety, alfy, muy, x, px, y, py, t, pt, dx, dpx, dy, dpy, wx, phix, dmux, wy, phiy, dmuy, ddx, ddpx, ddy, ddpy, r11, r12, r21, r22, energy, angle, k0l, k0sl, k1l, k1sl, k2l, k2sl, k3l, k3sl, k4l, k4sl, k5l, k5sl, k6l, k6sl, k7l, k7sl, k8l, k8sl, k9l, k9sl, k10l, k10sl, ksi, hkick, vkick, tilt, e1, e2, h1, h2, hgap, fint, fintx, volt, lag, freq, harmon, slot_id, assembly_id, mech_sep, kmax, kmin, calib, polarity, alfa, beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33, alfa11, alfa12, alfa13, alfa21, alfa22, disp1, disp2, disp3, disp4')
    madx_instance.twiss(sequence=sequence_name, file=file_out)

    return madx_instance.table.twiss.dframe()

########################################################################
# savebeta to save twiss at given location in given sequence
######################################################################## 
def cpymad_savebeta(madx_instance, cpymad_logfile, sequence_name, savebeta_label, location):
    log_string = '! cpymad_savebeta called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    madx_command = 'savebeta, label=' + str(savebeta_label) + ', place=' + str(location) + ', sequence=' + str(sequence_name) + ';'
    madx_instance.input(madx_command)
    print('cpymad_savebeta:: Remember to perform a Twiss to complete the savebeta')

########################################################################
# perform madx twiss with use sequence and beta0 return twiss dataframe
######################################################################## 
def cpymad_madx_twiss_beta0(madx_instance, cpymad_logfile, sequence_name, betazero, file_out=None):  
    if cpymad_check_and_use_sequence(madx_instance, cpymad_logfile, sequence_name):
        
        log_string = '! cpymad_madx_twiss_beta0 called for sequence ' + sequence_name + ' with betazero ' + betazero
        cpymad_write_to_logfile(cpymad_logfile, log_string)
        
        if file_out is None: file_out = sequence_name +'_madx_twiss.tfs'
    
        madx_instance.input('set, format="12.12f"')
        madx_instance.input('select, flag=twiss, column=keyword, name, s, l, betx, alfx, mux, bety, alfy, muy, x, px, y, py, t, pt, dx, dpx, dy, dpy, wx, phix, dmux, wy, phiy, dmuy, ddx, ddpx, ddy, ddpy, r11, r12, r21, r22, energy, angle, k0l, k0sl, k1l, k1sl, k2l, k2sl, k3l, k3sl, k4l, k4sl, k5l, k5sl, k6l, k6sl, k7l, k7sl, k8l, k8sl, k9l, k9sl, k10l, k10sl, ksi, hkick, vkick, tilt, e1, e2, h1, h2, hgap, fint, fintx, volt, lag, freq, harmon, slot_id, assembly_id, mech_sep, kmax, kmin, calib, polarity, alfa, beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33, alfa11, alfa12, alfa13, alfa21, alfa22, disp1, disp2, disp3, disp4')
        madx_instance.twiss(sequence=sequence_name,  beta0=betazero, file=file_out)
        
        return madx_instance.table.twiss.dframe()
    
########################################################################
# perform ptc twiss with use sequence, return twiss dataframe
########################################################################  
def cpymad_ptc_twiss(madx_instance, cpymad_logfile, sequence_name, file_out=None):  
    if cpymad_check_and_use_sequence(madx_instance, cpymad_logfile, sequence_name):
        
        log_string = '! cpymad_ptc_twiss called for sequence ' + sequence_name
        cpymad_write_to_logfile(cpymad_logfile, log_string)
    
        if file_out is None: file_out = sequence_name +'_ptc_twiss.tfs'
    
        madx_instance.input('ptc_create_universe')
        madx_instance.input('ptc_create_layout, time=false,model=2, method=6, nst=5, exact=true')
        madx_instance.input('set, format="12.12f"')
        madx_instance.input('ptc_twiss, closed_orbit, icase=56, no=4, slice_magnets')
        madx_instance.input('ptc_end')
        madx_instance.input('write, table=ptc_twiss, file='+file_out)
        
        ptc_twiss = tfs.read(file_out)
        pandas_dataframe_lowercase_columns(ptc_twiss)
        
        return ptc_twiss

########################################################################
# perform ptc twiss without use sequence, return twiss dataframe
########################################################################  
def cpymad_ptc_twiss_nocheck(madx_instance, cpymad_logfile, sequence_name, file_out=None):  
    if file_out is None: file_out = sequence_name +'_ptc_twiss.tfs'

    log_string = '! cpymad_ptc_twiss_nocheck called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)

    madx_instance.input('ptc_create_universe')
    madx_instance.input('ptc_create_layout, time=false,model=2, method=6, nst=5, exact=true')
    madx_instance.input('set, format="12.12f"')
    madx_instance.input('ptc_twiss, closed_orbit, icase=56, no=4, slice_magnets')
    madx_instance.input('ptc_end')
    madx_instance.input('write, table=ptc_twiss, file='+file_out)

    ptc_twiss = tfs.read(file_out)
    pandas_dataframe_lowercase_columns(ptc_twiss)
    
    return ptc_twiss
    
########################################################################
# perform madx plot to post-script file
########################################################################  
def cpymad_make_ps_plot(madx_instance, cpymad_logfile, sequence_name):
    log_string = '! cpymad_make_ps_plot called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
        
    madx_instance.input('setplot, font=4, xsize=34, ysize=25;')
    madx_instance.input('plot, table=twiss, haxis=s, vaxis1=betx,bety,dx, VAXIS2=mux,muy, vmin=0,vmax=30., interpolate, title='+str(sequence_name)+', style=100, file=plot;')
     
########################################################################
# perform seqedit to cycle starting position
########################################################################        
def cpymad_start_sequence_at(madx_instance, cpymad_logfile, sequence_name, cpymad_element):   
    log_string = '! cpymad_start_sequence_at called for sequence ' + sequence_name + ' start element ' + cpymad_element
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('flatten;')
    madx_instance.input(str('cycle, start='+cpymad_element+';'))
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')
    # return the name of the additionally created marker for easy removal later
    return str(sequence_name + cpymad_element + '_p_')
    
########################################################################
# perform seqedit to flatten (unpack) sequence
########################################################################
def cpymad_flatten_sequence(madx_instance, cpymad_logfile, sequence_name):
    log_string = '! cpymad_flatten_sequence called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)    
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('select, flag=seqedit, clear;')     
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')
     
########################################################################
# save sequence to file
########################################################################   
def cpymad_save_sequence(madx_instance, cpymad_logfile, sequence_name, savename, beam=False, bare=False, mad8=False):
    log_string = '! cpymad_save_sequence called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    madx_instance.input(str('save, sequence='+sequence_name+', file='+savename+', beam='+str(beam)+', bare='+str(bare)+', mad8='+str(mad8)+';'))   

########################################################################
# load sequence from file
########################################################################    
def cpymad_load_sequence_from_file(madx_instance, cpymad_logfile, sequence_file):
    log_string = '! cpymad_load_sequence_from_file called for file ' + sequence_file
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    # Find sequence in sequence file to store name as string
    f = open(sequence_file, "r")
    file_lines = f.readlines()
    seq_count = 0
    for line in file_lines:
        if ': sequence,' in line:
            print(line)
            seq_name = line.split(':')[0]
            if seq_count == 1: stored_name = line.split(':')[0]
            seq_count +=1
    if seq_count >= 2: print('cpymad_load_sequence_from_file::warning: multiple sequences in input file, first selected is ', stored_name)
    
    madx_instance.input(str('call, file='+str(sequence_file)+';'))
    sequence_name = cpymad_get_active_sequence_name(madx_instance)
    
    return sequence_name
    
########################################################################
# perform seqedit to insert element
########################################################################
def cpymad_insert_element(madx_instance, cpymad_logfile, sequence_name, cpymad_element, cpymad_class, cpymad_at, cpymad_from=None):
    log_string = '! cpymad_insert_element called for sequence ' + sequence_name + ' inserted element ' + cpymad_element
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    # Flatten lattice 
    # (unpack any subsequences within the sequence until the sequence is composed of a simple list of elements)
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('flatten;')
    if cpymad_from is None: madx_instance.input('INSTALL, ELEMENT='+str(cpymad_element)+', CLASS='+str(cpymad_class)+', AT='+str(cpymad_at)+';')
    else: madx_instance.input('INSTALL, ELEMENT='+str(cpymad_element)+', CLASS='+str(cpymad_class)+', AT='+str(cpymad_at)+', FROM='+str(cpymad_from)+';')
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')
    
########################################################################
# print element in sequence
########################################################################    
def cpymad_print_element_from_sequence(madx_instance, cpymad_logfile, cpymad_sequence, element):    
    if element in cpymad_sequence.elements:
        print(cpymad_sequence.elements[element])
        return True
    else:
        print('cpymad_print_element_from_sequence:: element',element,' not present in sequence ', cpymad_sequence.name)
        return False
        
########################################################################
# return element in sequence as string
########################################################################         
def cpymad_return_element_from_sequence(madx_instance, cpymad_logfile, cpymad_sequence, element):    
    if element in cpymad_sequence.elements: return str(cpymad_sequence.elements[element])    
    else: print('cpymad_print_element_from_sequence:: element',element,' not present in sequence ', cpymad_sequence.name)

########################################################################
# print global element
########################################################################    
def cpymad_print_global_element(madx_instance, cpymad_logfile, element):    
    if element in madx_instance.elements:
        print(madx_instance.elements[element])
        return True
    else:
        print('cpymad_print_global_element:: element',element,' not present in global MAD-X variables')
        return False
        
########################################################################
# return global element as string
########################################################################        
def cpymad_return_global_element(madx_instance, cpymad_logfile, element):    
    if element in madx_instance.elements: return str(madx_instance.elements[element])         
    else: print('cpymad_print_global_element:: element',element,' not present in global MAD-X variables')

########################################################################
# return detailed element string
########################################################################  
def cpymad_return_detailed_element(madx_instance, cpymad_logfile, element_name):   
    
    if element_name in madx_instance.elements: 
        #element_name = madx_instance.elements.element.name
        
        element_class = str(find_element_class(madx_instance.elements[element_name]))
        
        if element_class == ('sbend'):
            # store angle, k1 etc variables
            bend = str(element_name + ' : ' + element_class + ', l=' + str(madx_instance.elements[element_name].defs.l) + ', angle=' + str(madx_instance.elements[element_name].defs.angle) + ', k1=' + str(madx_instance.elements[element_name].defs.k1)+ ', e1=' + str(madx_instance.elements[element_name].defs.e1)+ ', e2=' + str(madx_instance.elements[element_name].defs.e2) + ', fint=' + str(madx_instance.elements[element_name].defs.fint) + ', fintx=' + str(madx_instance.elements[element_name].defs.fintx) + ';')
            #print(bend)
            return bend
        
        elif element_class == ('rbend'):
            # store angle, k1 etc variables
            bend = str(element_name + ' : ' + element_class + ', l=' + str(madx_instance.elements[element_name].defs.l) + ', angle=' + str(madx_instance.elements[element_name].defs.angle) + ', k1=' + str(madx_instance.elements[element_name].defs.k1)+ ', e1=' + str(madx_instance.elements[element_name].defs.e1)+ ', e2=' + str(madx_instance.elements[element_name].defs.e2) + ', fint=' + str(madx_instance.elements[element_name].defs.fint) + ', fintx=' + str(madx_instance.elements[element_name].defs.fintx) + ';')
            #print(bend)
            return bend
            
        elif element_class == ('quadrupole'):
            # store k1 variable
            quad = str(element_name + ' : ' + element_class + ', l=' + str(madx_instance.elements[element_name].defs.l) + ', k1=' + str(madx_instance.elements[element_name].defs.k1)+ ';')
            return quad
            
        elif element_class == ('sextupole'):
            # store k2 variable
            sext = str(element_name + ' : ' + element_class + ', l=' + str(madx_instance.elements[element_name].defs.l) + ', k2=' + str(madx_instance.elements[element_name].defs.k2)+ ';')
            return sext
            
        elif element_class == ('kicker'):
            kicker = str(element_name + ' : ' + element_class + ', kick=' + str(madx_instance.elements[element_name].defs.kick) + ';')
            return kicker
        
        elif element_class == ('hkicker'):
            kicker = str(element_name + ' : ' + element_class + ', kick=' + str(madx_instance.elements[element_name].defs.kick) + ';')
            return kicker
        
        elif element_class == ('vkicker'):
            kicker = str(element_name + ' : ' + element_class + ', kick=' + str(madx_instance.elements[element_name].defs.kick) + ';')
            return kicker
        
        elif element_class == ('drift'):
            drift = str(element_name + ' : ' + element_class + ', l=' + str(madx_instance.elements[element_name].defs.l) + ';')
            return drift
        
        elif element_class == ('marker'):
            drift = str(element_name + ' : ' + element_class + ';')
            return drift
            
        elif element_class == ('monitor'):
            monitor = str(element_name + ' : ' + element_class + ';')
            return monitor
            
        elif element_class == ('multipole'):
            multipole = str(madx_instance.elements[element_name])
            return multipole
            
        else:
            print('cpymad_return_detailed_element::warning: ', element_name, ' not found')
      
    else: 
        print('cpymad_print_global_element:: element',element_name,' not present in global MAD-X variables')

########################################################################
# simple duplicate and instantiate element (uses values not variables)
########################################################################  
def cpymad_create_element_duplicate(madx_instance, cpymad_logfile, old_element, new_element):
    log_string = '! cpymad_create_element_duplicate called for element ' + old_element
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    start_str = cpymad_return_global_element(madx_instance, cpymad_logfile, old_element)
    new_str = str(new_element) + ': ' + str(start_str).split(':')[1]
    madx_instance.input(new_str)
    return new_element
    
########################################################################
# complete duplicate and instantiate element (uses variables)
########################################################################  
def cpymad_create_complete_element_duplicate(madx_instance, cpymad_logfile, old_element, new_element):
    log_string = '! cpymad_create_complete_element_duplicate called for element ' + old_element
    cpymad_write_to_logfile(cpymad_logfile, log_string)     
    
    start_str = cpymad_return_detailed_element(madx_instance, cpymad_logfile, old_element)
    new_str = str(new_element) + ': ' + str(start_str).split(':')[1]
    print(new_str)
    madx_instance.input(new_str)
    return new_element

########################################################################
# select multiple elements for sequence editing
########################################################################  
def cpymad_select_multiple_elements(madx_instance, cpymad_logfile, sequence_name, element_list):
    log_string = '! cpymad_select_multiple_elements called for sequence ' + sequence_name + ', selected elements:'
    for el in element_list:
        log_string += '\n!' + el
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';')) 
    madx_instance.input('flatten;')
    madx_instance.input('select, flag=seqedit, clear;')    
    for el in element_list:        
        madx_instance.input(str('select, flag=seqedit, pattern=\"^'+el+'\";'))

########################################################################
# remove previously selected elements
########################################################################          
def cpymad_remove_selected_elements(madx_instance, cpymad_logfile, sequence_name):    
    log_string = '! cpymad_remove_selected_elements called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('flatten;')
    madx_instance.input('remove, element=selected')
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')
    
########################################################################
# replace previously selected elements
########################################################################   
def cpymad_replace_selected_elements(madx_instance, cpymad_logfile, sequence_name, new_element):    
    log_string = '! cpymad_replace_selected_elements called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';')) 
    madx_instance.input('flatten;')
    madx_instance.input(str('replace, element=selected, by='+str(new_element)+';'))
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')
    
########################################################################
# perform seqedit to replace element in a sequence
########################################################################      
def cpymad_replace_element_in_sequence(madx_instance, cpymad_logfile, sequence_name, old_element, new_element):
    
    log_string = '! cpymad_replace_element_in_sequence called for element ' +old_element+ ' in sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('flatten;')   
    madx_instance.input(str('replace, element='+str(old_element)+', by='+str(new_element)+';'))
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')

########################################################################
# perform seqedit to rename element in a sequence
########################################################################     
def cpymad_rename_element_in_sequence(madx_instance, cpymad_logfile, cpymad_sequence, old_element, new_element):
    
    log_string = '! cpymad_replace_element_in_sequence called for element ' +old_element+ ' in sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    # find the element
    if old_element in cpymad_sequence.elements:
        # create a copy of the element with a new name
        duplicate_element = cpymad_create_complete_element_duplicate(madx_instance, cpymad_logfile, old_element, new_element)
        
        # replace the element
        cpymad_replace_element_in_sequence(madx_instance, cpymad_logfile, cpymad_sequence.name, old_element, duplicate_element)
    
    else: print('cpymad_rename_element_in_sequence:: Requested element', old_element, ' not present in sequence ', cpymad_sequence.name)

########################################################################
# print all elements with the number of duplications in a given sequence 
########################################################################         
def cpymad_sequence_print_duplicate_names(madx_instance, cpymad_logfile, cpymad_sequence):
    element_counter = dict()
    for element in cpymad_sequence.elements:
        if element.name in element_counter.keys():
            element_counter[element.name] = element_counter[element.name] + 1
        else:
            element_counter[element.name] = 1
    
    for key, value in sorted(element_counter.items()):
        print(key, '->', value)

########################################################################
# perform a seqedit to extract a sequence - note cannot end with marker 
########################################################################
def cpymad_extract_sequence(madx_instance, cpymad_logfile, sequence_name, new_sequence_name, seq_start, seq_end):     
    log_string = '! cpymad_extract_sequence called for sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string) 
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('flatten;')   
    madx_instance.input(str('extract, sequence='+sequence_name+', from='+str(seq_start)+', to='+str(seq_end)+', newname='+new_sequence_name+';'))
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')
    
########################################################################
# perform a seqedit to remove an element from a sequence 
########################################################################
def cpymad_remove_element_from_sequence(madx_instance, cpymad_logfile, sequence_name, element):
    log_string = '! cpymad_remove_element_from_sequence called for element '+element+' in sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';'))
    madx_instance.input('flatten;')   
    madx_instance.input(str('remove, element='+element+';'))
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')  
    
########################################################################
# perform a seqedit to remove an element class from a sequence 
########################################################################
def cpymad_remove_all_elements_from_sequence(madx_instance, cpymad_logfile, sequence_name, element_class):      
    log_string = '! cpymad_remove_all_elements_from_sequence called for element class '+element_class+' in sequence ' + sequence_name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    madx_instance.input(str('seqedit, sequence='+sequence_name+';')) 
    madx_instance.input('select, flag=seqedit, clear;')      
    madx_instance.input('flatten;')       
    madx_instance.input(str('select, flag=seqedit, class='+element_class+';'))
    madx_instance.input('remove, element=selected;')
    madx_instance.input('flatten;')
    madx_instance.input('endedit;')  

########################################################################
# print all element names and the number of times duplicated in sequence
########################################################################
def cpymad_sequence_print_duplicate_names(madx_instance, cpymad_logfile, cpymad_sequence):
    element_counter = dict()
    for element in cpymad_sequence.elements:
        if element.name in element_counter.keys():
            element_counter[element.name] = element_counter[element.name] + 1
        else:
            element_counter[element.name] = 1
    
    for key, value in sorted(element_counter.items()):
        print(key, '->', value)

########################################################################
# return dict of all element names and duplication count in sequence
########################################################################
def cpymad_sequence_return_duplicate_names(madx_instance, cpymad_logfile, cpymad_sequence):
    element_counter = dict()
    for element in cpymad_sequence.elements:
        if element.name in element_counter.keys():
            element_counter[element.name] = element_counter[element.name] + 1
        else:
            element_counter[element.name] = 1
    
    return element_counter

########################################################################
# return dict of all element names and duplication count in sequence  
# not including drifts
########################################################################
def cpymad_sequence_return_duplicate_names_no_drifts(madx_instance, cpymad_logfile, cpymad_sequence):
    element_counter = dict()
    for element in cpymad_sequence.elements:
        if element.name in element_counter.keys():
            element_counter[element.name] = element_counter[element.name] + 1
        else:
            if ('drift' or 'DRIFT') not in str(element):
                element_counter[element.name] = 1
    
    return element_counter

########################################################################
# rename a single repeated element such that no names are duplicated 
########################################################################
def cpymad_individual_element_naming_in_sequence(madx_instance, cpymad_logfile, cpymad_sequence, repeated_element, new_name_start=None, verbose=True):
    log_string = '! cpymad_individual_element_naming_in_sequence called for element name '+repeated_element+' in sequence ' + cpymad_sequence.name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    count_iterator = 0
    
    if new_name_start is None: new_name_start = repeated_element
    
    for element in cpymad_sequence.elements:
        if verbose: print('Iterating through elements in ', cpymad_sequence.name, ' : ',element.name )
            
        if element.name == repeated_element:
            if verbose: print('Found element with name matching requirements: ', element)
                
            # copy, rename, instantiate a new instance in MAD-X
            new_element = cpymad_create_complete_element_duplicate(madx_instance, cpymad_logfile, repeated_element, str(new_name_start +'_'+str(count_iterator)))
            if verbose: print('Duplicated to: ', cpymad_return_detailed_element(madx_instance, cpymad_logfile, new_element))
            
            # replace
            cpymad_replace_element_in_sequence(madx_instance, cpymad_logfile, cpymad_sequence.name, repeated_element, new_element)
            print('Replacing element ', repeated_element, ' with ', new_element)
            count_iterator+=1
            
########################################################################
# rename all duplicates (not drifts) in sequence with appended _i to name 
########################################################################
def cpymad_rename_all_duplicates(madx_instance, cpymad_logfile, cpymad_sequence):
    log_string = '! cpymad_rename_all_duplicates called for sequence ' + cpymad_sequence.name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    duplicate_dict = cpymad_sequence_return_duplicate_names_no_drifts(madx_instance, cpymad_logfile, cpymad_sequence)
    
    for el_name, el_count in duplicate_dict.items():
    #for el_name in duplicate_dict.keys():
        if el_count >= 2:
            cpymad_individual_element_naming_in_sequence(madx_instance, cpymad_logfile, cpymad_sequence, el_name)
        
    final_dict = cpymad_sequence_return_duplicate_names_no_drifts(madx_instance, cpymad_logfile, cpymad_sequence)
    
    return final_dict

########################################################################
# rename all duplicates (including drifts) in sequence with appended _i to name 
########################################################################
def cpymad_rename_all_duplicates_inc_drifts(madx_instance, cpymad_logfile, cpymad_sequence):
    log_string = '! cpymad_rename_all_duplicates_inc_drifts called for sequence ' + cpymad_sequence.name
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    duplicate_dict = cpymad_sequence_return_duplicate_names(madx_instance, cpymad_logfile, cpymad_sequence)
    
    for el_name, el_count in duplicate_dict.items():
        if el_count >= 2:
            cpymad_individual_element_naming_in_sequence(madx_instance, cpymad_logfile, cpymad_sequence, el_name)
        
    final_dict = cpymad_sequence_return_duplicate_names(madx_instance, cpymad_logfile, cpymad_sequence)
    
    return final_dict
    
########################################################################
# find element class in cpymad description and return class string
########################################################################
def find_element_class(cpymad_element):    
    elements = ['quadrupole', 'sbend', 'marker', 'hkicker', 'vkicker', 'kicker', 'drift', 'marker', 'rbend', 'sextupole', 'octupole', 'multipole', 'solenoid', 'tkicker', 'rfcavity', 'rfmultipole', 'crabcavity', 'elseparator', 'monitor', 'hmonitor', 'vmonitor', 'instrument', 'placeholder', 'collimator', 'ecollimator', 'rcollimator', 'beambeam', 'matrix','yrotation', 'srotation','translation','changeref']
    cpymad_element = str(cpymad_element)    
    
    for el in elements:
        if el in cpymad_element:
            return(el)        
    
    print('find_element_class: No matching MAD-X class found in ', cpymad_element)    
    return False

########################################################################
# replace 2 elements with single element
########################################################################
def cpymad_replace_two_elements_with_element(madx_instance, cpymad_logfile, madx_sequence, element_list, new_madx_element):
    sequence_list, sequence_dict = cpymad_sequence_to_dict(madx_instance, cpymad_logfile, madx_sequence)
    sequence_list.append('fake_marker')
    
    log_string = '! cpymad_replace_two_elements_by_element called for sequence '+ str(madx_sequence.name) 

    # Iterate over list of element names
    for i in range(len(sequence_list)):
        if (element_list[0] in sequence_list[i]) and (element_list[1] in sequence_list[i+1]):
            print('cpymad_replace_two_elements_by_element: elements ' +sequence_list[i]+ ' and ' +sequence_list[i+1]+ ' replaced with ' + new_madx_element.name)
            
            # Insert new_element_name
            cpymad_insert_element(madx_instance, cpymad_logfile, madx_sequence.name, new_madx_element.name, find_element_class(new_madx_element), '0.0', sequence_list[i])
            cpymad_flatten_sequence(madx_instance, cpymad_logfile, madx_sequence.name)
            
            # Remove two elements
            cpymad_remove_element_from_sequence(madx_instance, cpymad_logfile, madx_sequence.name, sequence_list[i])
            cpymad_flatten_sequence(madx_instance, cpymad_logfile, madx_sequence.name)
            cpymad_remove_element_from_sequence(madx_instance, cpymad_logfile, madx_sequence.name, sequence_list[i+1])
            cpymad_flatten_sequence(madx_instance, cpymad_logfile, madx_sequence.name)
            
            log_string+='\n! elements ' +sequence_list[i]+ ' and ' +sequence_list[i+1]+ ' replaced with ' + new_madx_element.name

    cpymad_write_to_logfile(cpymad_logfile, log_string)

########################################################################
# Manual equivalent to saving sequence
########################################################################
def cpymad_save_detailed_line(madx_instance, cpymad_logfile, madx_sequence, save_name):
    log_string = '! cpymad_save_detailed_sequence called for sequence '+ str(madx_sequence.name)
    cpymad_write_to_logfile(cpymad_logfile, log_string)
    
    # Populate element dictionary and sequence list
    elements = dict()
    sequence = []
    for el in madx_sequence.elements:
        sequence.append(el.name)
        if el.name not in elements:
            elements[el.name] = cpymad_return_detailed_element(madx_instance, cpymad_logfile, el.name) 
            #print(elements[el.name])
    
    # Remove pesky None values from Dictionary
    filtered = {k: v for k, v in elements.items() if v is not None}
    elements.clear()
    elements.update(filtered)
    
    # Open file
    f = open(save_name, 'w')
    
    # First print madx globals (variables)
    ignore_globals = ['version', 'pi', 'twopi', 'degrad', 'raddeg', 'e', 'amu0', 'emass', 'mumass', 'nmass', 'pmass', 'clight', 'qelect', 'hbar', 'erad', 'prad', 'none', 'twiss_tol']
    for x, y in madx_instance.globals.defs.items():
        if x in ignore_globals:
            pass
        else:
            write_str =  str(x)+ ' = '+ str(y)+ ';\n'
            f.write(write_str)
        
    # Next print element dictionary
    for x, y in sorted(elements.items()):
        #print(y)
        write_str = '\n'+ str(y)
        f.write(write_str)
    
    # Sequence format requires AT positions
    #ring: sequence, l = 163.7413679;
    #d3_0, at = 0.125;
    #b1_0, at = 1.054614143;
    #f.write('\n', str(madx_sequence.name), ': sequence, l = ', str(madx_sequence.length) ,';')
    
    # Line format takes position based on element order in list
    #ring :line(sp_inj,3*sp)
    write_str = '\n'+str(madx_sequence.name)+': line=('
    f.write(write_str)
    
    # write sequence (note last element added separately to avoid comma)
    for x in sequence[:-1]: 
        write_str =str(x)+','
        f.write(write_str)
        
    f.write(sequence[-1])
    
    f.write(');')
    
    f.close()
########################################################################
########################################################################
#                         PLOTTING STUFF
########################################################################
########################################################################

########################################################################
# Block diagram
########################################################################
def block_diagram(ax1, df_myTwiss, ptc_twiss=False):
    
    # Remove Borders
    ax1.spines['top'].set_visible(False);
    ax1.spines['bottom'].set_visible(False);
    ax1.spines['left'].set_visible(False);
    ax1.spines['right'].set_visible(False);
    ax1.axes.get_yaxis().set_visible(False); 
    ax1.axes.get_xaxis().set_visible(False);  
    
    s_key =  's'
    keyword = 'keyword'  
    
    ###########
    ## QUADS ##
    ########### 
    if ptc_twiss: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-DF.iloc[i].l, 0.), DF.iloc[i].l, 1.0, color='r', alpha=0.5));
    
    ###############
    ## Sextupole ##
    ###############
    if ptc_twiss: key = 'SEXTUPOLE'
    else: key =  'sextupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-DF.iloc[i].l, 0.), DF.iloc[i].l, 1.0, color='green', alpha=0.5));

    ###########
    ## BENDS ##
    ########### 
    if ptc_twiss: key = 'SBEND'
    else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-DF.iloc[i].l, 0.), DF.iloc[i].l, 1.0, color='b', alpha=0.5));
            
    ############
    ## Kicker ##
    ############     
    kicker_length=0.5
    kicker_height = 1.0
    if ptc_twiss: key = 'KICKER'
    else: key =  'kicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        ax1.add_patch(patches.Rectangle( (DF.iloc[i].s, 0.), kicker_length, kicker_height, color='c', alpha=0.5));

    if ptc_twiss: key = 'HKICKER'
    else: key =  'hkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        ax1.add_patch(patches.Rectangle( (DF.iloc[i].s, 0.), kicker_length, kicker_height, color='c', alpha=0.5));

    if ptc_twiss: key = 'VKICKER'
    else: key =  'vkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        ax1.add_patch(patches.Rectangle( (DF.iloc[i].s, 0.), kicker_length, kicker_height, color='c', alpha=0.5));

    custom_lines = [Line2D([0], [0], color='b', lw=4, alpha=0.5),
                    Line2D([0], [0], color='r', lw=4, alpha=0.5),
                    Line2D([0], [0], color='green', lw=4, alpha=0.5),
                    Line2D([0], [0], color='cyan', lw=4, alpha=0.5)]

    if ptc_twiss:        
        ax1.set_xlim(0, df_myTwiss.headers['LENGTH']);
    else:
        ax1.set_xlim(0, df_myTwiss.iloc[-1].s);
    ax1.legend(custom_lines, ['Dipole', 'Quadrupole', 'Sextupole', 'Kicker']);

########################################################################
# Plot closed orbit
########################################################################
def cpymad_plot_CO(madx_instance, df_myTwiss, sequence_name, save_file, ptc_twiss=False, beta_rel=None):
        
    if ptc_twiss:
        if beta_rel == None:
            print('cpymad_plot_CO: Error: ptc_twiss selected but no beta_rel supplied for dispersion normalisation')
            return False
    
    # Plot title = sequence_name + tunes
    qx = madx_instance.table.summ.q1[0]
    qy = madx_instance.table.summ.q2[0]     
    plot_title = sequence_name +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
    
    # Start Plot
    heights = [1, 3, 2, 2]
    fig2 = plt.figure(figsize=(10,8),facecolor='w', edgecolor='k',constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig2, height_ratios=heights)
    
    # Block diagram
    f2_ax1 = fig2.add_subplot(spec2[0])
    f2_ax1.set_title(plot_title)  
    
    if ptc_twiss:
        block_diagram(f2_ax1, df_myTwiss, ptc_twiss=True)
    else: 
        block_diagram(f2_ax1, df_myTwiss, ptc_twiss=False)
    
    # Plot betas   
    f2_ax2 = fig2.add_subplot(spec2[1], sharex=f2_ax1)  
    f2_ax2.plot(df_myTwiss['s'], df_myTwiss['betx'],'b', label='$\\beta_x$')
    f2_ax2.plot(df_myTwiss['s'], df_myTwiss['bety'],'r', label='$\\beta_y$')    
    
    f2_ax2.legend(loc=2)
    f2_ax2.set_ylabel(r'$\beta_{x,y}$[m]')
    f2_ax2.grid(which='both', ls=':', lw=0.5, color='k')
    #f2_ax2.set_xlabel('s [m]')
    #f2_ax2.set_xticklabels([])
    
    if np.min(df_myTwiss['bety']) < np.min(df_myTwiss['betx']): bet_min = round_down_n(np.min(df_myTwiss['bety']),5)
    else: bet_min = round_down_n(np.min(df_myTwiss['betx']),5)
    if np.max(df_myTwiss['bety']) > np.max(df_myTwiss['betx']): bet_max = round_up_n(np.max(df_myTwiss['bety']),10)
    else: bet_max = round_up_n(np.max(df_myTwiss['betx']),10)        
    f2_ax2.set_ylim(bet_min,bet_max)
    
    ax2 = f2_ax2.twinx()   # instantiate a second axes that shares the same x-axis
    if ptc_twiss:     
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp1']/beta_rel,'green', label='$D_x$')
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp3']/beta_rel,'purple', label='$D_y$')
        key_dx = 'disp1';        key_dy = 'disp3';  
    
    else:
        ax2.plot(df_myTwiss['s'], df_myTwiss['dx'],'green', label='$D_x$')
        ax2.plot(df_myTwiss['s'], df_myTwiss['dy'],'purple', label='$D_y$')
        key_dx = 'dx';        key_dy = 'dy';  
        
    ax2.legend(loc=1)
    ax2.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down_n(np.min(df_myTwiss[key_dy]),1)
    else: d_min = round_down_n(np.min(df_myTwiss[key_dx]),1)    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_n(np.max(df_myTwiss[key_dy]),10)
    else: d_max = round_up_n(np.max(df_myTwiss[key_dx]),10) 
    ax2.set_ylim(d_min,d_max)
   
    
    f2_ax3 = fig2.add_subplot(spec2[2], sharex=f2_ax1)
    f2_ax3.plot(df_myTwiss['s'], df_myTwiss['x']*1E3,'k', lw=1.5, label='$x$')      
    #f2_ax3.legend(loc=2)
    f2_ax3.set_ylabel('x [mm]')
    f2_ax3.grid(which='both', ls=':', lw=0.5, color='k')
    
    f2_ax4 = fig2.add_subplot(spec2[3], sharex=f2_ax1)
    f2_ax4.plot(df_myTwiss['s'], df_myTwiss['y']*1E3,'k', lw=1.5, label='$y$')   
    f2_ax4.set_ylabel('y [mm]')
    f2_ax4.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss['y']) < np.min(df_myTwiss['x']): co_min = round_down_n(np.min(df_myTwiss['y']),10)
    else: co_min = round_down_n(np.min(df_myTwiss['x']),10)
    if np.max(df_myTwiss['y']) > np.max(df_myTwiss['x']): co_max = round_up_n(np.max(df_myTwiss['y']),10)
    else: co_max = round_up_n(np.max(df_myTwiss['x']),10)   
        
    f2_ax4.set_ylim(co_min,co_max)
    f2_ax3.set_ylim(co_min,co_max)
    
    f2_ax4.set_xlabel('s [m]')
    
    #f2_ax4 = fig2.add_subplot(spec2[4], sharex=f2_ax1)   
    if save_file != None: plt.savefig(save_file)

########################################################################
# Plot closed orbit correction
########################################################################
def cpymad_plot_CO_correction(madx_instance, df_myTwiss, df_myTwiss2, sequence_name, save_file, ptc_twiss=False, beta_rel=None):
        
    if ptc_twiss:
        if beta_rel == None:
            print('cpymad_plot_CO: Error: ptc_twiss selected but no beta_rel supplied for dispersion normalisation')
            return False
    
    # Plot title = sequence_name + tunes
    qx = madx_instance.table.summ.q1[0]
    qy = madx_instance.table.summ.q2[0]     
    plot_title = sequence_name +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
    
    # Start Plot
    heights = [1, 3, 2, 2]
    fig2 = plt.figure(figsize=(10,8),facecolor='w', edgecolor='k',constrained_layout=True);
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig2, height_ratios=heights);
    
    # Block diagram
    f2_ax1 = fig2.add_subplot(spec2[0]);
    f2_ax1.set_title(plot_title);  
    
    if ptc_twiss:
        block_diagram(f2_ax1, df_myTwiss, ptc_twiss=True)
    else: 
        block_diagram(f2_ax1, df_myTwiss, ptc_twiss=False)
    
    # Plot betas   
    f2_ax2 = fig2.add_subplot(spec2[1], sharex=f2_ax1);  
    f2_ax2.plot(df_myTwiss['s'], df_myTwiss['betx'],'b', label='$\\beta_x$');
    f2_ax2.plot(df_myTwiss['s'], df_myTwiss['bety'],'r', label='$\\beta_y$');  
    
    f2_ax2.legend(loc=2);
    f2_ax2.set_ylabel(r'$\beta_{x,y}$[m]');
    f2_ax2.grid(which='both', ls=':', lw=0.5, color='k');
    #f2_ax2.set_xlabel('s [m]')
    #f2_ax2.set_xticklabels([])
    
    if np.min(df_myTwiss['bety']) < np.min(df_myTwiss['betx']): bet_min = round_down_n(np.min(df_myTwiss['bety']),5)
    else: bet_min = round_down_n(np.min(df_myTwiss['betx']),5)
    if np.max(df_myTwiss['bety']) > np.max(df_myTwiss['betx']): bet_max = round_up_n(np.max(df_myTwiss['bety']),10)
    else: bet_max = round_up_n(np.max(df_myTwiss['betx']),10)        
    f2_ax2.set_ylim(bet_min,bet_max);
    
    ax2 = f2_ax2.twinx();   # instantiate a second axes that shares the same x-axis
    if ptc_twiss:     
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp1']/beta_rel,'green', label='$D_x$');
        ax2.plot(df_myTwiss['s'], df_myTwiss['disp3']/beta_rel,'purple', label='$D_y$');
        key_dx = 'disp1';        key_dy = 'disp3';  
    
    else:
        ax2.plot(df_myTwiss['s'], df_myTwiss['dx'],'green', label='$D_x$');
        ax2.plot(df_myTwiss['s'], df_myTwiss['dy'],'purple', label='$D_y$');
        key_dx = 'dx';        key_dy = 'dy';  
        
    ax2.legend(loc=1);
    ax2.set_ylabel(r'$D_{x,y}$ [m]', color='green');  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor='green');
    ax2.grid(which='both', ls=':', lw=0.5, color='green');

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down_n(np.min(df_myTwiss[key_dy]),1)
    else: d_min = round_down_n(np.min(df_myTwiss[key_dx]),1)    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_n(np.max(df_myTwiss[key_dy]),10)
    else: d_max = round_up_n(np.max(df_myTwiss[key_dx]),10) 
    ax2.set_ylim(d_min,d_max);
   
    # Include standard deviations (RMS in mm) of orbits on plots 
    if np.std(df_myTwiss2['x']) == 0.0:
        print('cpymad_plot_CO_correction: Error: Uncorrected closed orbit is not perturbed, please apply errors')
        return False
    std_x = helper_functions.round_sig(np.std(df_myTwiss2['x'])*1E3,3)
    std_x_corr = helper_functions.round_sig(np.std(df_myTwiss['x'])*1E3,3)
    txt_x = r'$x_{RMS}$ = ' + str(std_x) + ' mm\n' + r'$x_{corrected~RMS}$ = '+ str(std_x_corr)+ ' mm'
    std_y = helper_functions.round_sig(np.std(df_myTwiss2['y'])*1E3,3)
    std_y_corr = helper_functions.round_sig(np.std(df_myTwiss['y'])*1E3,3)
    txt_y = r'$y_{RMS}$ = ' + str(std_y) + ' mm\n' + r'$y_{corrected~RMS}$ = '+ str(std_y_corr)+ ' mm'
    
    f2_ax3 = fig2.add_subplot(spec2[2], sharex=f2_ax1);
    f2_ax3.plot(df_myTwiss2['s'], df_myTwiss2['x']*1E3,'r', lw=1.5, label='$x$');
    f2_ax3.plot(df_myTwiss['s'], df_myTwiss['x']*1E3,'k', lw=1.5, label='$x_{corrected}$');     
    #f2_ax3.legend(loc=2)
    f2_ax3.set_ylabel('x [mm]');
    f2_ax3.grid(which='both', ls=':', lw=0.5, color='k');
    f2_ax3.text(np.max(df_myTwiss['s'])/2, 5, txt_x, fontsize=10);
    
    f2_ax4 = fig2.add_subplot(spec2[3], sharex=f2_ax1);
    f2_ax4.plot(df_myTwiss2['s'], df_myTwiss2['y']*1E3,'r', lw=1.5, label='$y$')  ;
    f2_ax4.plot(df_myTwiss['s'], df_myTwiss['y']*1E3,'k', lw=1.5, label='$y_{corrected}$');   
    f2_ax4.set_ylabel('y [mm]');
    f2_ax4.grid(which='both', ls=':', lw=0.5, color='k');
    f2_ax4.text(np.max(df_myTwiss['s'])/2, 5, txt_y, fontsize=10);
    
    if np.min(df_myTwiss['y']) < np.min(df_myTwiss['x']): co_min = round_down_n(np.min(df_myTwiss['y']),10)
    else: co_min = round_down_n(np.min(df_myTwiss['x']),10)
    if np.max(df_myTwiss['y']) > np.max(df_myTwiss['x']): co_max = round_up_n(np.max(df_myTwiss['y']),10)
    else: co_max = round_up_n(np.max(df_myTwiss['x']),10)   
    
    f2_ax3.legend(loc=2);
    f2_ax4.legend(loc=2);
    
    f2_ax3.set_ylim(co_min,co_max);
    f2_ax4.set_ylim(co_min,co_max);
    
    f2_ax4.set_xlabel('s [m]');
    
    #f2_ax4 = fig2.add_subplot(spec2[4], sharex=f2_ax1)   
    if save_file != None: 
        plt.savefig(save_file);
        plt.close();       
        
########################################################################
# auxilliary functions for plotting function 
########################################################################
def plotLatticeSeries(ax,series, height=1., v_offset=0., color='r',alpha=0.5):
    aux=series
    ax.add_patch(
    patches.Rectangle(
        (aux.s-aux.l, v_offset-height/2.),   # (x,y)
        aux.l,          # width
        height,          # height
        color=color, alpha=alpha
    )
    )
    return;

def plotLatticeSeriesCaps(ax,series, height=1., v_offset=0., color='r',alpha=0.5):
    aux=series
    ax.add_patch(
    patches.Rectangle(
        (aux.S-aux.L, v_offset-height/2.),   # (x,y)
        aux.L,          # width
        height,          # height
        color=color, alpha=alpha
    )
    )
    return;

def round_up_10(x):
    return int(math.ceil(x / 10.0)) * 10

def round_down_10(x):
    return int(math.floor(x / 10.0)) * 10

def round_up_p1(x):
    return int(math.ceil(x / 0.1)) * 0.1

def round_down_p1(x):
    return int(math.floor(x / 0.1)) * 0.1

def round_down(x):
    return int(math.floor(x / 1.0)) * 1

def round_up_n(x,n):
    return int(math.ceil(x / 10.0)) * int(n)

def round_down_n(x,n):
    return int(math.floor(x / 10.0)) * int(n)


########################################################################
# automatically plot madx twiss
########################################################################
def cpymad_plot_madx_twiss(madx_instance, df_myTwiss, title=None, savename=None, use_caps=False):
        
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    
    if use_caps: 
        qx = madx_instance.table.summ.Q1[0]
        qy = madx_instance.table.summ.Q2[0]
    else:
        qx = madx_instance.table.summ.q1[0]
        qy = madx_instance.table.summ.q2[0]               
    
    if title is None:        
        active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
        plot_title = active_seq +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
    else:
        plot_title = title + ' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
        

    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)      
    plt.title(plot_title) 
    
    if use_caps is False:
        try:    plt.plot(df_myTwiss['s'], 0*df_myTwiss['s'],'k')
        except KeyError: 
            try: 
                plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')
                use_caps = True
                print('cpymad_plotTwiss::use_caps = True')
            except: print('cpymad_plotTwiss::unkown bug')
    else:
        plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')        

    if use_caps: 
        s_key = 'S'
        keyword = 'KEYWORD'        
    else: 
        s_key =  's'
        keyword = 'keyword'   
    
    quad_max = 0.
    dipole_max = 0.
    
    if use_caps: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            plotLatticeSeriesCaps(plt.gca(),aux, height=aux.K1L, v_offset=aux.K1L/2, color='r')
            if np.max(abs(aux.K1L)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.K1L)))
        else: 
            plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
            if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))
    
    if use_caps: key = 'MULTIPOLE' 
    else: key =  'multipole' 
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            plotLatticeSeriesCaps(plt.gca(),aux, height=aux.K1L, v_offset=aux.K1L/2, color='r')
            if np.max(abs(aux.K1L)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.K1L)))
        else: 
            plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
            if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))


    #plt.ylim(-.065,0.065)
    color = 'red'
    ax1.set_ylabel('1/f=K1L [m$^{-1}$]', color=color)  # we already handled the x-label with ax1
    ax1.tick_params(axis='y', labelcolor=color)
    plt.grid()
    
    plt.ylim(-quad_max,quad_max)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis   
    
    color = 'blue'
    ax2.set_ylabel('$\\theta$=K0L [rad]', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)

    if use_caps: key = 'SBEND'
    else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            plotLatticeSeriesCaps(plt.gca(),aux, height=aux.ANGLE, v_offset=aux.ANGLE/2, color='b')
            if np.max(abs(aux.ANGLE)) > dipole_max: dipole_max = round_up_p1(np.max(abs(aux.ANGLE)))
        else: 
            plotLatticeSeries(plt.gca(),aux, height=aux.angle, v_offset=aux.angle/2, color='b')
            if np.max(abs(aux.angle)) > dipole_max: dipole_max = round_up_p1(np.max(abs(aux.angle)))

    plt.ylim(-dipole_max,dipole_max)

    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    if use_caps: key_betx = 'BETX';        key_bety = 'BETY';
    else:        key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel('s [m]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    if use_caps: key_dx = 'DX';        key_dy = 'DY';
    else:        key_dx = 'dx';        key_dy = 'dy';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx],'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy],'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    #plt.ylim(round_down(np.min(df_myTwiss[key_dx])), round_up_10(np.max(df_myTwiss[key_dx])))
     
    if savename is None: pass
    else: plt.savefig(savename)


########################################################################
# automatically plot ptc twiss
########################################################################
def cpymad_plot_ptc_twiss(madx_instance, df_myTwiss, title=None, savename=None, use_caps=True):
    
    pandas_dataframe_lowercase_columns(df_myTwiss)
    
    gamma_key = 'GAMMA'; pc_key='PC';
    ptc_twiss_read_Header = dict(df_myTwiss.headers)
    gamma_rel = ptc_twiss_read_Header[gamma_key]
    beta_rel = np.sqrt( 1. - (1./gamma_rel**2) )
    p_mass_GeV = 0.93827208816 #Proton mass GeV
    tot_energy = gamma_rel * p_mass_GeV
    kin_energy = tot_energy - p_mass_GeV
    momentum = ptc_twiss_read_Header[pc_key]

    print('Relativistic Gamma = ', round(gamma_rel,3))
    print('Relativistic Beta = ', round(beta_rel,3))
    print('Total Energy = ', round(tot_energy,4), 'GeV')
    print('Kinetic Energy = ', round(kin_energy*1E3,3), 'MeV')
    print('momentum = ', round(momentum,3), 'GeV/c')
    
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)     
    
    plt.plot(df_myTwiss['s'], 0*df_myTwiss['s'],'k')
        
    s_key =  's'
    keyword = 'keyword'   
    Q1 = ptc_twiss_read_Header['Q1']
    Q2 = ptc_twiss_read_Header['Q2']
      
    if title is None:
        active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
        plot_title = active_seq +' Q1='+format(Q1,'2.3f')+', Q2='+ format(Q2,'2.3f')
    else:
        plot_title = title + ' Q1='+format(Q1,'2.3f')+', Q2='+ format(Q2,'2.3f')
    
    plt.title(plot_title)
            
    quad_max = 0.
    dipole_max = 0.
    
    if use_caps: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]

        plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
        if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))
    
    if use_caps: key = 'MULTIPOLE' 
    else: key =  'multipole' 
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
 
        plotLatticeSeries(plt.gca(),aux, height=aux.k1l, v_offset=aux.k1l/2, color='r')
        if np.max(abs(aux.k1l)) > quad_max: quad_max = round_up_p1(np.max(abs(aux.k1l)))

    color = 'red'
    ax1.set_ylabel('1/f=K1L [m$^{-1}$]', color=color)  # we already handled the x-label with ax1
    ax1.tick_params(axis='y', labelcolor=color)
    plt.grid()
    plt.ylim(-quad_max,quad_max)
    
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'blue'
    ax2.set_ylabel('$\\theta$=K0L [rad]', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)

    if use_caps: key = 'SBEND'
    else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]

        plotLatticeSeries(plt.gca(),aux, height=aux.angle, v_offset=aux.angle/2, color='b')
        if np.max(abs(aux.angle)) > dipole_max: dipole_max = round_up_p1(np.max(abs(aux.angle)))

    plt.ylim(-dipole_max,dipole_max)
    
    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel('s [m]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    key_dx = 'disp1';        key_dy = 'disp3';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx]/beta_rel,'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy]/beta_rel,'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')
    
    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    
    if savename is None: pass
    else: plt.savefig(savename)
    
########################################################################
# automatically plot madx twiss block with s as x-axis
########################################################################
def cpymad_plot_madx_twiss_block(madx_instance, df_myTwiss, title=None, savename=None, use_caps=False):
        
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    
    if use_caps: 
        qx = madx_instance.table.summ.Q1[0]
        qy = madx_instance.table.summ.Q2[0]
    else:
        qx = madx_instance.table.summ.q1[0]
        qy = madx_instance.table.summ.q2[0]               
    
    if title is None:        
        active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
        plot_title = active_seq +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
    else:
        plot_title = title + ' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
        
    # Start Plot
    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)      
    plt.title(plot_title) 
    
    # Empty plot
    if use_caps is False:
        try:    plt.plot(df_myTwiss['s'], 0*df_myTwiss['s'],'k')
        except KeyError: 
            try: 
                plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')
                use_caps = True
                print('cpymad_plotTwiss::use_caps = True')
            except: print('cpymad_plotTwiss::unkown bug')
    else:
        plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')        

    if use_caps: 
        s_key = 'S'
        keyword = 'KEYWORD'        
    else: 
        s_key =  's'
        keyword = 'keyword'   
    
    #---------START BLOCK DIAGRAM-------------
    
    ###########
    ## QUADS ##
    ###########    
    if use_caps: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-DF.iloc[i].L, 0.), DF.iloc[i].L, 1.0, color='r', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-DF.iloc[i].l, 0.), DF.iloc[i].l, 1.0, color='r', alpha=0.5))
    
    ###############
    ## Sextupole ##
    ###############
    if use_caps: key = 'SEXTUPOLE'
    else: key =  'sextupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-DF.iloc[i].L, 0.), DF.iloc[i].L, 1.0, color='green', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-DF.iloc[i].l, 0.), DF.iloc[i].l, 1.0, color='green', alpha=0.5))

    ###########
    ## BENDS ##
    ########### 
    if use_caps: key = 'SBEND'
    else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-DF.iloc[i].L, 0.), DF.iloc[i].L, 1.0, color='b', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-DF.iloc[i].l, 0.), DF.iloc[i].l, 1.0, color='b', alpha=0.5))
            
    ############
    ## Kicker ##
    ############     
    kicker_length=0.5
    kicker_height = 1.0
    
    if use_caps: key = 'KICKER'
    else: key =  'kicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    if use_caps: key = 'HKICKER'
    else: key =  'hkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    if use_caps: key = 'VKICKER'
    else: key =  'vkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    custom_lines = [Line2D([0], [0], color='b', lw=4, alpha=0.5),
                    Line2D([0], [0], color='r', lw=4, alpha=0.5),
                    Line2D([0], [0], color='green', lw=4, alpha=0.5),
                    Line2D([0], [0], color='cyan', lw=4, alpha=0.5)]

    ax1.legend(custom_lines, ['Dipole', 'Quadrupole', 'Sextupole', 'Kicker'])
            
    
    #---------START TWISS-------------    
    
    ###########
    ## TWISS ##
    ###########

    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    if use_caps: key_betx = 'BETX';        key_bety = 'BETY';
    else:        key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel('s [m]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    if use_caps: key_dx = 'DX';        key_dy = 'DY';
    else:        key_dx = 'dx';        key_dy = 'dy';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx],'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy],'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    #plt.ylim(round_down(np.min(df_myTwiss[key_dx])), round_up_10(np.max(df_myTwiss[key_dx])))
     
    if savename is None: pass
    else: plt.savefig(savename)
    
########################################################################
# automatically plot madx twiss block with x phase as x-axis
########################################################################
def cpymad_plot_madx_twiss_block_phase_x(madx_instance, df_myTwiss, title=None, savename=None, use_caps=False):
        
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    
    if use_caps: 
        qx = madx_instance.table.summ.Q1[0]
        qy = madx_instance.table.summ.Q2[0]
    else:
        qx = madx_instance.table.summ.q1[0]
        qy = madx_instance.table.summ.q2[0]               
    
    if title is None:        
        active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
        plot_title = active_seq +' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
    else:
        plot_title = title + ' Q1='+format(qx,'2.3f')+', Q2='+ format(qy,'2.3f')
        
    # Start Plot
    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)      
    plt.title(plot_title) 
    
    # Empty plot
    if use_caps is False:
        try:    plt.plot(df_myTwiss['mux'], 0*df_myTwiss['mux'],'k')
        except KeyError: 
            try: 
                plt.plot(df_myTwiss['MUX'], 0*df_myTwiss['MUX'],'k')
                use_caps = True
                print('cpymad_plotTwiss::use_caps = True')
            except: print('cpymad_plotTwiss::unkown bug')
    else:
        plt.plot(df_myTwiss['S'], 0*df_myTwiss['S'],'k')        

    if use_caps: 
        s_key = 'MUX'
        keyword = 'KEYWORD'        
    else: 
        s_key =  'mux'
        keyword = 'keyword'   
    
    #---------START BLOCK DIAGRAM-------------
    kicker_length= .01
    kicker_height = 1.0    
    
    ###########
    ## QUADS ##
    ###########    
    if use_caps: key = 'QUADRUPOLE'
    else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MUX, 0.), kicker_length, 1.0, color='r', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mux, 0.), kicker_length, 1.0, color='r', alpha=0.5))
    
    ###############
    ## Sextupole ##
    ###############
    if use_caps: key = 'SEXTUPOLE'
    else: key =  'sextupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MUX, 0.), kicker_length, 1.0, color='green', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mux, 0.), kicker_length, 1.0, color='green', alpha=0.5))

    ###########
    ## BENDS ##
    ########### 
    if use_caps: key = 'SBEND'
    else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MUX, 0.), kicker_length, 1.0, color='b', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mux, 0.), kicker_length, 1.0, color='b', alpha=0.5))
            
    ############
    ## Kicker ##
    ############     
    if use_caps: key = 'KICKER'
    else: key =  'kicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MUX, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mux, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    if use_caps: key = 'HKICKER'
    else: key =  'hkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MUX, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mux, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    if use_caps: key = 'VKICKER'
    else: key =  'vkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MUX, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mux, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    custom_lines = [Line2D([0], [0], color='b', lw=4, alpha=0.5),
                    Line2D([0], [0], color='r', lw=4, alpha=0.5),
                    Line2D([0], [0], color='green', lw=4, alpha=0.5),
                    Line2D([0], [0], color='cyan', lw=4, alpha=0.5)]

    ax1.legend(custom_lines, ['Dipole', 'Quadrupole', 'Sextupole', 'Kicker'])
    
    #---------START TWISS-------------    
    
    ###########
    ## TWISS ##
    ###########

    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    if use_caps: key_betx = 'BETX';        key_bety = 'BETY';
    else:        key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel(r'$\mu_x$ [2$\pi$]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    if use_caps: key_dx = 'DX';        key_dy = 'DY';
    else:        key_dx = 'dx';        key_dy = 'dy';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx],'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy],'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    #plt.ylim(round_down(np.min(df_myTwiss[key_dx])), round_up_10(np.max(df_myTwiss[key_dx])))
     
    if savename is None: pass
    else: plt.savefig(savename)               
    
########################################################################
# automatically plot ptc twiss block diagram with s as x-axis
########################################################################
def cpymad_plot_ptc_twiss_block(madx_instance, df_myTwiss, title=None, savename=None, use_caps=False):
    
    pandas_dataframe_lowercase_columns(df_myTwiss)
    
    gamma_key = 'GAMMA'; pc_key='PC';
    ptc_twiss_read_Header = dict(df_myTwiss.headers)
    gamma_rel = ptc_twiss_read_Header[gamma_key]
    beta_rel = np.sqrt( 1. - (1./gamma_rel**2) )
    p_mass_GeV = 0.93827208816 #Proton mass GeV
    tot_energy = gamma_rel * p_mass_GeV
    kin_energy = tot_energy - p_mass_GeV
    momentum = ptc_twiss_read_Header[pc_key]

    print('Relativistic Gamma = ', round(gamma_rel,3))
    print('Relativistic Beta = ', round(beta_rel,3))
    print('Total Energy = ', round(tot_energy,4), 'GeV')
    print('Kinetic Energy = ', round(kin_energy*1E3,3), 'MeV')
    print('momentum = ', round(momentum,3), 'GeV/c')
    
    # Start Plot
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)     
    
    plt.plot(df_myTwiss['s'], 0*df_myTwiss['s'],'k')
    
    s_key =  's'
    keyword = 'keyword'   
    Q1 = ptc_twiss_read_Header['Q1']
    Q2 = ptc_twiss_read_Header['Q2']
      
    if title is None:
        active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
        plot_title = active_seq +' Q1='+format(Q1,'2.3f')+', Q2='+ format(Q2,'2.3f')
    else:
        plot_title = title + ' Q1='+format(Q1,'2.3f')+', Q2='+ format(Q2,'2.3f')
    
    plt.title(plot_title)    
    
    #---------START BLOCK DIAGRAM-------------
    element_length = 0.05
    
    ###########
    ## QUADS ##
    ###########      
    #if use_caps: 
    key = 'QUADRUPOLE'
    #else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-element_length, 0.), element_length, 1.0, color='r', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-element_length, 0.), element_length, 1.0, color='r', alpha=0.5))
    
    ###############
    ## Sextupole ##
    ###############
    #if use_caps: 
    key = 'SEXTUPOLE'
    #else: key =  'sextupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-element_length, 0.), element_length, 1.0, color='green', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-element_length, 0.), element_length, 1.0, color='green', alpha=0.5))

    ###########
    ## BENDS ##
    ########### 
    #if use_caps: 
    key = 'SBEND'
    #else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-element_length, 0.), element_length, 1.0, color='b', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-element_length, 0.), element_length, 1.0, color='b', alpha=0.5))
            
    ############
    ## Kicker ##
    ############     
    kicker_length=0.5
    kicker_height = 1.0
    
    #if use_caps: 
    key = 'KICKER'
    #else: key =  'kicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-element_length, 0.), element_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-element_length, 0.), element_length, kicker_height, color='c', alpha=0.5))

    #if use_caps: 
    key = 'HKICKER'
    #else: key =  'hkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-element_length, 0.), element_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-element_length, 0.), element_length, kicker_height, color='c', alpha=0.5))

    #if use_caps: 
    key = 'VKICKER'
    #else: key =  'vkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].S-element_length, 0.), element_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].s-element_length, 0.), element_length, kicker_height, color='c', alpha=0.5))

    custom_lines = [Line2D([0], [0], color='b', lw=4, alpha=0.5),
                    Line2D([0], [0], color='r', lw=4, alpha=0.5),
                    Line2D([0], [0], color='green', lw=4, alpha=0.5),
                    Line2D([0], [0], color='cyan', lw=4, alpha=0.5)]

    ax1.legend(custom_lines, ['Dipole', 'Quadrupole', 'Sextupole', 'Kicker'])
    
    #---------START TWISS-------------    
    
    ###########
    ## TWISS ##
    ###########

    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel('s [m]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    key_dx = 'disp1';        key_dy = 'disp3';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx],'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy],'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    #plt.ylim(round_down(np.min(df_myTwiss[key_dx])), round_up_10(np.max(df_myTwiss[key_dx])))
     
    if savename is None: pass
    else: plt.savefig(savename)             
    
########################################################################
# automatically plot ptc twiss block diagram with x phase as x-axis
########################################################################
def cpymad_plot_ptc_twiss_block_phase_x(madx_instance, df_myTwiss, title=None, savename=None, use_caps=False):
    
    pandas_dataframe_lowercase_columns(df_myTwiss)
    
    gamma_key = 'GAMMA'; pc_key='PC';
    ptc_twiss_read_Header = dict(df_myTwiss.headers)
    gamma_rel = ptc_twiss_read_Header[gamma_key]
    beta_rel = np.sqrt( 1. - (1./gamma_rel**2) )
    p_mass_GeV = 0.93827208816 #Proton mass GeV
    tot_energy = gamma_rel * p_mass_GeV
    kin_energy = tot_energy - p_mass_GeV
    momentum = ptc_twiss_read_Header[pc_key]

    print('Relativistic Gamma = ', round(gamma_rel,3))
    print('Relativistic Beta = ', round(beta_rel,3))
    print('Total Energy = ', round(tot_energy,4), 'GeV')
    print('Kinetic Energy = ', round(kin_energy*1E3,3), 'MeV')
    print('momentum = ', round(momentum,3), 'GeV/c')
    
    # Start Plot
    fig = plt.figure(figsize=(13,8),facecolor='w', edgecolor='k')
    ax1=plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)     
    
    plt.plot(df_myTwiss['mu1'], 0*df_myTwiss['mu1'],'k')
    
    s_key =  'mu1'
    keyword = 'keyword'   
    Q1 = ptc_twiss_read_Header['Q1']
    Q2 = ptc_twiss_read_Header['Q2']
      
    if title is None:
        active_seq = str(cpymad_get_active_sequence(madx_instance)).split('\'')[1]
        plot_title = active_seq +' Q1='+format(Q1,'2.3f')+', Q2='+ format(Q2,'2.3f')
    else:
        plot_title = title + ' Q1='+format(Q1,'2.3f')+', Q2='+ format(Q2,'2.3f')
    
    plt.title(plot_title)
    
    #---------START BLOCK DIAGRAM-------------
    kicker_length= .005
    kicker_height = 1.0    
    
    ###########
    ## QUADS ##
    ###########    
    #if use_caps: 
    key = 'QUADRUPOLE'
    #else: key =  'quadrupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MU1, 0.), kicker_length, 1.0, color='r', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mu1, 0.), kicker_length, 1.0, color='r', alpha=0.5))
    
    ###############
    ## Sextupole ##
    ###############
    #if use_caps: 
    key = 'SEXTUPOLE'
    #else: key =  'sextupole'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MU1, 0.), kicker_length, 1.0, color='green', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mu1, 0.), kicker_length, 1.0, color='green', alpha=0.5))

    ###########
    ## BENDS ##
    ########### 
    #if use_caps: 
    key = 'SBEND'
    #else: key =  'sbend'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        if use_caps: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MU1, 0.), kicker_length, 1.0, color='b', alpha=0.5))
        else: 
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mu1, 0.), kicker_length, 1.0, color='b', alpha=0.5))
            
    ############
    ## Kicker ##
    ############     
    #if use_caps: 
    key = 'KICKER'
    #else: key =  'kicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MU1, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mu1, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    #if use_caps: 
    key = 'HKICKER'
    #else: key =  'hkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MU1, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mu1, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    #if use_caps: 
    key = 'VKICKER'
    #else: key =  'vkicker'
    DF=df_myTwiss[(df_myTwiss[keyword]==key)]
    for i in range(len(DF)):
        aux=DF.iloc[i]
        
        if use_caps:              
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].MU1, 0.), kicker_length, kicker_height, color='c', alpha=0.5))
        else:  
            ax1.add_patch(patches.Rectangle( (DF.iloc[i].mu1, 0.), kicker_length, kicker_height, color='c', alpha=0.5))

    custom_lines = [Line2D([0], [0], color='b', lw=4, alpha=0.5),
                    Line2D([0], [0], color='r', lw=4, alpha=0.5),
                    Line2D([0], [0], color='green', lw=4, alpha=0.5),
                    Line2D([0], [0], color='cyan', lw=4, alpha=0.5)]

    ax1.legend(custom_lines, ['Dipole', 'Quadrupole', 'Sextupole', 'Kicker'])
    
    #---------START TWISS-------------    
    
    ###########
    ## TWISS ##
    ###########

    # large subplot
    plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2,sharex=ax1)
    key_betx = 'betx';        key_bety = 'bety';        
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_betx],'b', label='$\\beta_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_bety],'r', label='$\\beta_y$')
    plt.legend(loc=2)
    plt.ylabel(r'$\beta_{x,y}$[m]')
    plt.xlabel(r'$\mu_x$ [2$\pi$]')
    plt.grid(which='both', ls=':', lw=0.5, color='k')
    
    if np.min(df_myTwiss[key_bety]) < np.min(df_myTwiss[key_betx]): bet_min = round_down_10(np.min(df_myTwiss[key_bety]))
    else: bet_min = round_down_10(np.min(df_myTwiss[key_betx]))
    if np.max(df_myTwiss[key_bety]) > np.max(df_myTwiss[key_betx]): bet_max = round_up_10(np.max(df_myTwiss[key_bety]))
    else: bet_max = round_up_10(np.max(df_myTwiss[key_betx]))        
    plt.ylim(bet_min,bet_max)

    ax3 = plt.gca().twinx()   # instantiate a second axes that shares the same x-axis
    key_dx = 'disp1';        key_dy = 'disp3';  
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dx],'green', label='$D_x$')
    plt.plot(df_myTwiss[s_key], df_myTwiss[key_dy],'purple', label='$D_y$')
    ax3.legend(loc=1)
    ax3.set_ylabel(r'$D_{x,y}$ [m]', color='green')  # we already handled the x-label with ax1
    ax3.tick_params(axis='y', labelcolor='green')
    plt.grid(which='both', ls=':', lw=0.5, color='green')

    if np.min(df_myTwiss[key_dy]) < np.min(df_myTwiss[key_dx]): d_min = round_down(np.min(df_myTwiss[key_dy]))
    else: d_min = round_down(np.min(df_myTwiss[key_dx]))    
    if np.max(df_myTwiss[key_dy]) > np.max(df_myTwiss[key_dx]): d_max = round_up_10(np.max(df_myTwiss[key_dy]))
    else: d_max = round_up_10(np.max(df_myTwiss[key_dx]))        
    plt.ylim(d_min,d_max)
    #plt.ylim(round_down(np.min(df_myTwiss[key_dx])), round_up_10(np.max(df_myTwiss[key_dx])))
     
    if savename is None: pass
    else: plt.savefig(savename)               
