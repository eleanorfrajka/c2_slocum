import os
import datetime

def create_figdir():
    #figdir = '../results/figfiles/yyyy-mm/';
    today = datetime.datetime.now()
    yyyymm = str(today.year)+'-'+str(today.month)
    figdir = os.path.join('..','03-results','figures',yyyymm)
    if not os.path.isdir(figdir):
        os.makedirs(figdir, exist_ok=True)
        #os.mkdir(figdir)
    return figdir
        
def save_figure(fig,figname):
    figdir = create_figdir()
    figname_with_ext = 'fig_'+figname+'.png'
    fig.savefig(os.path.join(figdir,figname_with_ext))

def cat_data_path(data_folder,filename):
    input_file = os.path.join('..',data_folder,filename)
    return input_file

def cat_raw_path(fname):
    output_file_with_path = os.path.join('..','01-data','01-raw',fname)
    return output_file_with_path

def cat_interim_path(fname):
    output_file_with_path = os.path.join('..','01-data','02-intermediate',fname)
    return output_file_with_path

def cat_proc_path(fname):
    output_file_with_path = os.path.join('..','01-data','03-processed',fname)
    return output_file_with_path