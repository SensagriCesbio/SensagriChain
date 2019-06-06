#!/usr/bin/python
# Ludo 06/06/2019
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import matplotlib.backends.backend_pdf as mpdf

def PlotCM(matList,OAList,IntList,RecList,PreList,FSList,datesList,classnames,filepdf,title,method):
    # Colormap
    cmap = 'rainbow'
    filename = "%s-%s.pdf"%(filepdf,method)
    pdfplot = mpdf.PdfPages(filename)
 
    for k,mat in enumerate(matList):
	# Adapt values shape
	date = datesList[k]
        OA = OAList[k]
        Int = IntList[k]
	Rec = RecList[k]
	Pre = PreList[k]
	FS = FSList[k] 

        Nmat = len(mat)  
        Rec = Rec.reshape(1,-1).T
        Pre = [Pre]
        FS = FS.reshape(1,-1).T

        # Initialize plot
	plottitle = "%s - Average Confusion Matrix - %s - %s"%(title,method,date)
        fig = plt.figure(figsize=(20,15))
        fig.suptitle(plottitle, size = 24) 
        grid = plt.GridSpec(2, 3,width_ratios=[3, 1, 1], height_ratios=[3, 1])
        
        # Add axes
        ax1 = fig.add_subplot(grid[0,0])
        ax2 = fig.add_subplot(grid[0,1])
        ax3 = fig.add_subplot(grid[1,0])
        ax4 = fig.add_subplot(grid[0,2])
        
        # Add titles
        ax1.set_title("Overall Accuracy = %.2f +/- %.2f %%"%(OA,Int),size=20)
        ax2.set_title("Recall (%)")
        ax3.set_title("Precision (%)")
        ax4.set_title("FScore (%)")
        
        lin = np.arange(Nmat)
        linsup = np.arange(Nmat+1) + 0.5
        
        # Major ticks
        ax1.set_xticks(lin);
        ax1.set_yticks(lin);
        ax2.set_yticks(lin);
        ax3.set_xticks(lin);
        ax4.set_yticks(lin);
        
        # Minor ticks
        ax1.set_xticks(linsup, minor=True);
        ax1.set_yticks(linsup, minor=True);
        ax2.set_yticks(linsup, minor=True);
        ax3.set_xticks(linsup, minor=True);
        ax4.set_yticks(linsup, minor=True);
        
        # Labels for major ticks
        ax1.set_xticklabels([]);
        ax1.set_yticklabels([]);
        ax2.set_xticklabels([]);
        ax2.set_yticklabels([]);
        ax3.set_xticklabels([]);
        ax3.set_yticklabels([]);
        ax4.set_xticklabels([]);
        ax4.set_yticklabels([]);
        
        # Labels for minor ticks
        ax1.set_xticklabels(classnames, rotation=90, minor = True);
        ax1.set_yticklabels(classnames, minor = True);
        
        ax2.set_yticklabels(classnames, minor = True);
        ax3.set_xticklabels(classnames, rotation=90, minor = True);
        ax4.set_yticklabels(classnames, minor = True);
        
        # Gridlines based on minor ticks
        ax1.grid(which='major', linestyle = '-', color = 'black', linewidth=1.5)
        ax2.grid(which='major', linestyle = '-', axis = 'y', color = 'black')
        ax3.grid(which='major', linestyle = '-', axis = 'x', color= 'black')
        ax4.grid(which='major', linestyle = '-', axis = 'y', color = 'black')
        
        # col and row plot size
        extentMat = [0,len(mat[0]),len(mat[0]),0]
        extentCol = [0,1,len(mat[0]),0]
        extentRow = [0,len(mat[0]),1,0]
        
        for i in range(Nmat):
          for j in range(Nmat):
            #ax1.text(i+0.5,j+0.5,"%.1E"%(mat[i,j]), color='black', ha='center', va='center', rotation=45, size = 5)
            val = mat[j,i]
            if val < 10000:
              fval = "%d"%(val)
            else:
              e = int(np.floor(np.log10(abs(val))))
              m = val/10**e 
              fval = "%.1fe%d"%(m,e)
        
            ax1.text(i+0.5,j+0.5,fval, color='black', ha='center', va='center', rotation=45, size = 11)
        pl1 = ax1.imshow(mat, extent = extentMat, interpolation='none', aspect='equal',cmap = cmap, norm=LogNorm(vmin=0.1, vmax=5000000))
        
        for i in range(Nmat):
          for j in range(Nmat):
            #ax1.text(i+0.5,j+0.5,"%.1E"%(mat[i,j]), color='black', ha='center', va='center', rotation=45, size = 5)
            ax2.text(0.5,j+0.5,"%d"%(Rec[j][0]), color='black', ha='center', va='center', size = 10)
            ax3.text(i+0.5,0.5,"%d"%(Pre[0][i]), color='black', ha='center', va='center', size = 10)
            ax4.text(0.5,j+0.5,"%d"%(FS[j][0]), color='black', ha='center', va='center', size = 10)
        
        
        ax2.imshow(Rec, extent = extentCol, interpolation='none', aspect='equal', cmap = cmap)
        ax3.imshow(Pre, extent = extentRow, interpolation='none', aspect='equal', cmap = cmap)
        ax4.imshow(FS,  extent = extentCol, interpolation='none', aspect='equal', cmap = cmap)
        
        fig.colorbar(pl1, ax = ax1)
        #pdfplot.savefig(fig,bbox_inches='tight')
        pdfplot.savefig(fig)
    pdfplot.close()   

    #plt.tight_layout()
