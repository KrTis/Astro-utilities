import matplotlib.pyplot as plt
from math import ceil
import numpy as np
from uncertainties import ufloat
from collections import OrderedDict
from copy import deepcopy
from Tkinter import *
import os
import shutil
import ttk
import tkFileDialog
import tkMessageBox
from astroML_problematic import binned_statistic_2d, sigmaG, convert_to_stdev
from json import dump as jdump
from scipy.optimize import curve_fit
from scipy.stats import norm, moment
############################################################################

class comp:
	@classmethod
	def reportwin(self,direktorij,namel,report,reportl):
	    save_report=open(str(direktorij+'/report.tex'),'w')
	    save_report.write(report)
	    save_report.close()
	    plt.figure(figsize=(8,6))
	    ax=plt.gca()
	    plt.axis('off')
	    plt.table(cellText=reportl, colLabels=namel,loc='center')
	    plt.savefig(str(direktorij+'/report.png'))
	    plt.close()

class App:
	'''
	the Main program
	'''
	sizes=10
	defaulthist=50
	defaulthist2d=100
	defaultlabel=30
	labele0=OrderedDict([('$\\chi^2$', 1), ('$T_{per}$', 3), ('$e$', 4), ('$\\omega$', 5), ('$K_1$', 6), ('$K_2$', 7)])
	cmap_names=sorted(['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap',
		 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r',
		 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 
		'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r',
		 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 
		'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 
		'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary',
		 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper',
		 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 
		'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 
		'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 
		'hsv_r', 'jet', 'jet_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'prism', 'prism_r', 
		'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', 'summer_r', 
		'terrain', 'terrain_r', 'winter', 'winter_r'],key=str.lower)
	
	def __init__(self,master):
		self.Number_of_rows=0
		self.Number_of_columns=0
		self.chosen_cmap = StringVar(value='jet')
		self.f=''
		self.histbins=StringVar(value=App.defaulthist)
		self.histbins2d=StringVar(value=App.defaulthist2d)
		self.label_size=StringVar(value=App.defaultlabel)
		self.output_folder=StringVar(value='results')
		self.txtv=[StringVar() for j in range(self.Number_of_columns)]
		self.txtl=[StringVar() for j in range(self.Number_of_columns)]
		self.txtu=[StringVar() for j in range(self.Number_of_columns)]
		self.addv=StringVar()
		self.plot_output=IntVar(value=0)

  
		self.menubar = Menu(master)
		self.filemenu = Menu(self.menubar, tearoff=0)
		self.filemenu.add_command(label="Open",command=self.openf)
		self.filemenu.add_command(label="Output folder",command=self.o_folder)
		self.filemenu.add_separator()
		self.filemenu.add_command(label='Open limits', command=self.open_limits)
		self.filemenu.add_command(label='Save limits', command=self.save_limits)
		self.filemenu.add_command(label='Properties', command=self.properties_window)
		self.filemenu.add_separator()
		self.filemenu.add_command(label='Exit', command=lambda:master.quit())
		
				
	
		self.helpmenu=Menu(self.menubar, tearoff=0)
		self.helpmenu.add_command(label="About",command=self.clickAbout)
		self.helpmenu.add_command(label="How to use?",command=self.Help)		
	
		self.menubar.add_cascade(label="File", menu=self.filemenu)
		self.menubar.add_command(label="Run",command=self.start)
		self.menubar.add_cascade(label="Help", menu=self.helpmenu)
		master.config(menu=self.menubar)

	def openf(self):
		dialog_output=tkFileDialog.askopenfilename(parent=master,title='Choose a file')
		self.f=open(dialog_output)
		self.Number_of_rows=sum([1 for sumacija in self.f])
		self.f.seek(0)
		self.Number_of_columns=max([len(sumacija.strip().split()) for sumacija in self.f])
		self.f.seek(0)
		self.txtv=[StringVar() for j in range(self.Number_of_columns)]
		self.txtl=[StringVar() for j in range(self.Number_of_columns)]
		self.txtu=[StringVar() for j in range(self.Number_of_columns)]
		auxil=[[] for j in range(self.Number_of_columns)]
		self.f.seek(0)
		for line in self.f:
			line=line.split()
			for j in range(len(auxil)):
				auxil[j].append(float(line[j]))
		for j in range(self.Number_of_columns):
			self.txtl[j].set(min(auxil[j]))
			self.txtu[j].set(max(auxil[j]))

		separator = Frame(master,bd=1, relief=SUNKEN)
		separator.grid(row=1,column=0)
		additional=Frame(master,bd=1,relief=SUNKEN)
		additional.grid(row=2,column=0)
		Label(additional,text='Start offset').grid(row=0,column=0)
	
		Entry(additional,width=App.sizes,textvariable=self.addv).grid(row=0,column=1)

		Label(separator ,text='n').grid(row=1,column=0)
		Label(separator ,text='Names').grid(row=1,column=1)
		Label(separator ,text='Lower cuts').grid(row=1,column=2)
		Label(separator ,text='Upper cuts').grid(row=1,column=3)
		for k,j in App.labele0.items():
				if(j<self.Number_of_columns):		
					self.txtv[j].set(k)
		for i in range(self.Number_of_columns):
			Label(separator ,text=str(i)).grid(row=2+i,column=0)
			trt=0
			b=Entry(separator ,width=App.sizes,textvariable=self.txtv[i])
			b.grid(row=2+i,column=1)
			b=Entry(separator ,width=App.sizes,textvariable=self.txtl[i])
			b.grid(row=2+i,column=2)
			b=Entry(separator,width=App.sizes,textvariable=self.txtu[i])
			b.grid(row=2+i,column=3)

	def open_limits(self):	
		save_name=tkFileDialog.askopenfilename(parent=master,title='Open limits file')
		filesave=open(save_name,'r')

		j=0
		try:
			for line in filesave:
					line1=line.strip()
					line2=line1.split(':')
					self.txtv[j].set(line2[1])
					self.txtl[j].set(line2[2])
					self.txtu[j].set(line2[3])
					j=j+1
		except IndexError:
			if len(self.txtv)==0:
				tkMessageBox.showerror(message="No data file loaded, cannot import limits!")
			else:
				tkMessageBox.showerror(message="Limits not applicable!")
            		return
		filesave.close()
	
	def save_limits(self):
		save_name=tkFileDialog.asksaveasfilename(parent=master,title='Save limits file')
		filesave=open(save_name,'w')
		for i in range(self.Number_of_columns):
			filesave.write(str(str(i)+':'+str(self.txtv[i].get())+':'+str(self.txtl[i].get())+':'+str(self.txtu[i].get())+'\n'))
		filesave.close()
	
	def o_folder(self):	
		of=tkFileDialog.askdirectory(parent=master,title='Choose output folder')
		self.output_folder.set(of)	
		if not os.path.exists(of):
		   		 os.makedirs(of)
	def clickAbout(self):
	    toplevel=Toplevel(master)
	    toplevel.wm_title("About")
	    ABOUT_TEXT="Multi Dimensional Analyzer (2/6/2016)\n Made using Numpy, PyPlot and AstroML. \n Please report bugs to kresimir.tisanic@gmail.com"
	    about = Label(toplevel, text=ABOUT_TEXT, height=0, width=50)
	    about.pack()
	
	def Help(self):
	    toplevel=Toplevel(master)
	    toplevel.wm_title("How To Use?")
	    ABOUT_TEXT="MDA Takes an ASCII file as input.\n Columns should be separated by tabs or spaces. \n Upon hitting Run, it computes the statistic \n on parameters defined as columns."
	    about = Label(toplevel, text=ABOUT_TEXT, height=0, width=50)
	    about.pack()
	def properties_window(self):
		toplevel=Toplevel(master)
		toplevel.wm_title('Properties')
		self.colorframe=ttk.LabelFrame(toplevel,text='Colorbar')
		self.colorframe.grid(row=0,column=0)
		opt_cmap=ttk.OptionMenu(self.colorframe,self.chosen_cmap , 'jet', *App.cmap_names)
		opt_cmap.grid(row=0,column=0)
		Button(self.colorframe,text='Preview',command=self.show_colorbar).grid(row=0,column=1)
		self.plotframe=ttk.LabelFrame(toplevel,text='Plotting Options')
		self.plotframe.grid(row=1,column=0)
		Label(self.plotframe,text='Histogram bins').grid(row=0,column=0)
		Label(self.plotframe,text='2D Histogram bins').grid(row=1,column=0)
		Label(self.plotframe,text='Label sizes').grid(row=2,column=0)
		Entry(self.plotframe,width=App.sizes,textvariable=self.histbins).grid(row=0,column=1)
		Entry(self.plotframe,width=App.sizes,textvariable=self.histbins2d).grid(row=1,column=1)
		Entry(self.plotframe,width=App.sizes,textvariable=self.label_size).grid(row=2,column=1)
		Checkbutton(self.plotframe, text='Save plotted data', variable=self.plot_output).grid(row=3,column=0)
	def show_colorbar(self):
     
		fig=plt.figure(figsize=(30,30))

		sizep=int(ceil(np.sqrt(len(App.cmap_names))))
		ax = plt.subplot2grid((sizep, sizep),(0,0),colspan=sizep, rowspan=1)

		gradient = np.linspace(0, 1, 256)
		gradient = np.vstack((gradient, gradient))
		ax.set_title('Chosen colormap: '+self.chosen_cmap.get())
		ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(self.chosen_cmap.get()))
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)

		for x in range(1,len(App.cmap_names)+1):
		                  ax=fig.add_subplot(1+sizep, sizep,x+sizep)
		                  ax.set_title(App.cmap_names[x-1])
		                  ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(App.cmap_names[x-1]))
		                  ax.get_xaxis().set_visible(False)
		                  ax.get_yaxis().set_visible(False)
		fig.subplots_adjust(hspace=1,wspace=0.1,top=None, bottom=None, left=None, right=None)
		plt.show(fig)
	def indexed(self,label):
                             out={}
                             c=0
                             for x in label.keys():
                                 out[x]=c
                                 c=c+1
                             return out, c
	def gauss(self,x,a,b,c):
         return a*np.exp((x-b)**2/(2*c))
	def start(self):
		plt.ioff()
		#prerequisites
		if not self.f:
			return
		if not os.path.exists(self.output_folder.get()):
	   		 os.makedirs(self.output_folder.get())
		### definitions
		finishing_string='Name & Mean & Median & $\Sigma" & "K"\\\\ \n'
		starting_list=['Name','Mean','Median','Skewness', 'Kurtosis']
		finishing_list=[]
		histbins=int(self.histbins.get())
		histbins2d=int(self.histbins2d.get())
		labelsiz=int(self.label_size.get())
		textval=[a.get() for a in self.txtv]
		textu=[a.get() for a in self.txtu]
		textl=[a.get() for a in self.txtl]
		labele=OrderedDict([(textval[j],j) for j in range(len(textval)) if(len(textval[j])!=0)])
		labeleu=OrderedDict([(textval[j],float(textu[j])) for j in range(len(textval)) if((len(textu[j])!=0) and (len(textval[j])!=0))])
		labelel=OrderedDict([(textval[j],float(textl[j])) for j in range(len(textval)) if((len(textl[j])!=0) and (len(textval[j])!=0))])
		try:
			addit=int(self.addv.get())
		except ValueError:
			addit=0

		self.cmap_multicolor=plt.get_cmap(self.chosen_cmap.get())
		### reading data
		self.f.seek(0)
		X=np.array([[0. for i in range(self.Number_of_columns)] for j in range(self.Number_of_rows)])
		i=0
		j=0
		for line in self.f:
				if j>addit:
					line0=line.strip()
					line1=line0.split()
					passing=1
					for dicts, vals in labele.items():
						if dicts in labelel:
							if(labelel[dicts]>float(line1[vals])):
								passing=0
						if dicts in labeleu:
							if(labeleu[dicts]<float(line1[vals])):
								passing=0
					if(passing==1):			
						for j in range(self.Number_of_columns):
							X[i][j]=float(line1[j])
				
						i=i+1
				j=j+1
		X=X[0:i,:]
		X=np.array([[X[k][l] for l in range(self.Number_of_columns)] for k in range(i) ])
		if self.plot_output.get()==1:
                    with open(str(self.output_folder.get()+'/data.txt'),'w') as output_file:
                        np.savetxt(output_file,X)
                    
		
                        
		### Plotting histograms
		fig=plt.figure(figsize=(25,19))
		count=1
		dimension=len(labele)
		for ime, val in labele.items():
				ax=fig.add_subplot(2,int(np.ceil(dimension/2)),count)
				ax.set_xlabel(ime, size=labelsiz)
				ax.hist(X[:,val], histbins, normed=True, histtype='stepfilled',color='blue',facecolor='blue')
				param = norm.fit(X[:,val])
				moments= [moment(X[:,val], moment=i) for i in range(2,5)]
				sigmam=moments[0]**(0.5)
				skewness=moments[1]/moments[0]**(3/2.)
				kurtosis=moments[2]/moments[0]**(2.)-3.
				x=np.linspace(min(X[:,val]),max(X[:,val]),100)
				pdf_fitted = norm.pdf(x, *param)
				ax.plot(x, pdf_fitted, color='r')
				ax.axvline(np.mean(X[:,val]), color='orange', linestyle='--')
				ax.axvline(np.median(X[:,val]), color='green', linestyle='--')
				ax.xaxis.major.formatter._useMathText = True
				ax.ticklabel_format(style='sci', axis='x', scilimits=(-5,5))
				text='Mean and median:\n'+\
    "mean$\\rightarrow$ $ {:.2uL}$\n".format(ufloat(np.mean(X[:,val]),np.std(X[:,val])))+\
    "median$\\rightarrow$ $ {:.2uL}$".format(ufloat(np.median(X[:,val]),sigmaG(X[:,val])))+"\n"+\
    'Gaussian fit\n'+"peak $\\rightarrow\, {:.2uL}$".format(ufloat(param[0],param[1]))+"\n"+\
    "Moments:\n"+"$\Sigma=$ {:0.1f}\n".format(skewness)+\
    "$K=$ {:0.1f}".format(kurtosis)
				finishing_string=finishing_string+ime+" & "+"$ {:.2uL}$".format(ufloat(np.mean(X[:,val]),np.std(X[:,val])))+" & $ {:.2uL}$".format(ufloat(np.median(X[:,val]),sigmaG(X[:,val])))+"& $ {:0.1f}$".format(skewness)+"&"+"$ {:0.1f}$".format(kurtosis)+"\\\\ \n"			
				finishing_list.append([ime,"$ {:.2uL}$".format(ufloat(np.mean(X[:,val]),np.std(X[:,val])))," $ {:.2uL}$".format(ufloat(np.median(X[:,val]),sigmaG(X[:,val]))),"$ {:0.1f}$".format(skewness),"$ {:0.1f}$".format(kurtosis)])
				ax.text(.05,.8,text,transform = ax.transAxes,size='small')
				count=count+1
		plt.tight_layout()
		plt.savefig(self.output_folder.get()+'/Histograms.png')
		plt.close()
		comp.reportwin(self.output_folder.get(),starting_list,finishing_string,finishing_list)
  
		### Specific plots
		for names0,vals0 in labele.items():
				fig=plt.figure(figsize=(30,30))
				labele1=deepcopy(labele)

				del labele1[names0]
				labele2=deepcopy(labele1)

				counts=1
				indices1,dimension1=self.indexed(labele1)
				
				for names,vals in labele1.items():
					del labele2[names]
					indices2,dimension2=self.indexed(labele2)
					for names1, vals1 in labele2.items():
 		                   
						N0, xedges0, yedges0 = binned_statistic_2d(X[:,vals], X[:,vals1], X[:,labele[names0]], 'mean', bins=histbins2d)
						ax = plt.subplot2grid((dimension1,dimension1-1),(indices1[names1],indices1[names]),colspan=1, rowspan=1)
						im=ax.imshow(N0.T, origin='lower',extent=[xedges0[0], xedges0[-1], yedges0[0], yedges0[-1]], aspect='auto', interpolation='nearest', cmap=self.cmap_multicolor)
						plt.xlim(xedges0[0], xedges0[-1])
						plt.ylim(yedges0[0], yedges0[-1])
						if     indices1[names1]<dimension1-1:
           						ax.set_xticklabels('')
						else:
           						plt.xlabel(names, size=labelsiz)
           						ax.xaxis.major.formatter._useMathText = True
           						ax.ticklabel_format(style='sci', axis='x', scilimits=(-5,5))
						if     indices1[names]>0:
           						ax.set_yticklabels('')
						else:
           						plt.ylabel(names1, size=30)
           						ax.yaxis.major.formatter._useMathText = True
           						ax.ticklabel_format(style='sci', axis='y', scilimits=(-5,5))

						
						ax.yaxis.major.formatter._useMathText = True
						
                                  
						H, xbins, ybins = np.histogram2d(X[:,vals], X[:,vals1],bins=100)

						Nsigma = convert_to_stdev(np.log(H))
						
						cont=plt.contour(0.5 * (xbins[1:] + xbins[:-1]),0.5 * (ybins[1:] + ybins[:-1]),Nsigma.T,levels=[0.6827,0.6827,0.9545, 0.9545], colors=['.25','.25','0.5','0.5'],linewidths=2)					
						
						counts=counts+1
						
				
				self.cmap_multicolor.set_bad('w', 1.)
				fig.subplots_adjust(bottom=0.12,hspace=0,wspace=0)
				cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.025])
				cb=fig.colorbar(im, cax=cbar_ax, format=r'$%.1f$',orientation='horizontal')
				cb.set_label(str('$\\langle '+names0.replace('$','')+'\\rangle $'), size=labelsiz)
				
				plt.savefig(self.output_folder.get()+'/'+''.join([i for i in names0 if (i.isalpha() or i.isdigit())])+'.png',bbox_inches='tight')
				plt.close()




master = Tk()
App(master)

master.wm_title("Multi Dimensional Analyzer")
master.mainloop()



