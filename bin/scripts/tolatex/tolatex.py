#!/usr/bin/python
# -*- coding: utf-8 -*-
# tolatex: to construct latex file with python
# Ludo 29/06/2017 #
import os
import sys
import subprocess

class latexdoc:
    """
       Output to latex format
    """
    # Constructor #
    def __init__(self,outputfile,verbose=False):
        os.system("mkdir -p latex.output")
	directory = os.environ["SENSAGRICHAIN_HOME"] +"/bin/scripts/tolatex/" 
        self.outputfile = outputfile + ".tex"
	self.verbose = verbose
        with open(directory + "template.tex","r") as f:
            self.template = f.readlines()
            f.close()
        self.f = open(self.outputfile,"w")
       
        for i in range(len(self.template) - 1): 
            self.f.write(self.template[i])

    # Destructor #
    def __del__(self):
        if not self.f.closed:
            self.close()                

    # write #    
    def write(self,string):
        self.f.write(string + "\n")

    # title     
    def title(self,title,author = ""):
         self.f.write("\\title{" + title +"}\n")
         if len(author) > 0:
             self.f.write("\\author{%s}"%(author))
         self.f.write("\\date{\oldstylenums{\\today}}") 
         self.f.write("\maketitle\n")

    # tableau 
    def tableau(self,data,tformat,multicol = [],size = [], caption = ""):
       # if( (len(Cdata)!=len(TabData[0][0])+1) or (len(Rdata)!=datasize) ):
       #     print ""
       #     print "    Erreur in todata.tableau()"
       #     print "    Input data do not have matching dimenssions."
       #     print ""
       #     quit()
        self.write("\\begin{table}[h]")
        self.write("%\\scalebox{0.75}{")
        self.write("\\begin{tabular}{%s}"%(tformat))
        
        if len(multicol) > 0:
            #print "MULTI",len(multicol)
            self.write("\\hline")
            s = ""
            for m in range(len(multicol)):
                if m == 0:
                   s1 = ""; s2 = "|"; s3 = "|"
                elif m == len(multicol) - 1:
                   s1 = " & "; s2 = ""; s3 = "|"
                else:
                   s1 = " & "; s2 = ""; s3 = "||"
                s = s + "%s\\multicolumn{%d}{%sc%s}{%s}"%(s1,size[m],s2,s3,multicol[m])
            self.write(s)
            self.write("\\\\")
            self.write("\\hline")
        #else:
            #print "NO MULTICOL"

	first = False
        cidx = 0
        for mat in data:
            if first:
                self.write("\\hline")
            first = True
            self.write("\\hline")
            for l in mat:
                if type(l[0]) == str:
                    ptype = "%s "
                elif isinstance(l[0],int):
                    ptype = "%d "
                else:
                    ptype = "%.2f "
                if cidx == 0:
                    cidx = 1
                    colorcmd = "\\rowcolor[gray]{1}"
                else:
                    cidx = 0
                    colorcmd = "\\rowcolor[gray]{0.8}"
                s = colorcmd + ptype%(l[0])
                #self.write("\\hline")

                for c in range(1,len(l)):
                    if type(l[c]) == str:
                        ptype = "& %s "
                    elif isinstance(l[c],int):
                        ptype = "& %d "
                    else:
                        ptype = "& %.2f "
                    s = s + ptype%(l[c])
                s = s + " \\\\"
                self.write(s)
                      
        self.write("""\\hline
\\end{tabular}
%%}
\\caption{%s}
\\end{table}
"""%(caption))

    # matrix  
    def matrix(self,data,tformat,multicol = [],size = [], caption = ""):
       # if( (len(Cdata)!=len(TabData[0][0])+1) or (len(Rdata)!=datasize) ):
       #     print ""
       #     print "    Error in todata.tableau()"
       #     print "    Input data do not have matching dimenssions."
       #     print ""
       #     quit()
        self.write("\\begin{table}[H]")
        self.write("\\begin{tabular}{%s}"%(tformat))
        
        if len(multicol) > 0:
            #print "MULTI",len(multicol)
            self.write("\\hline")
            s = ""
            for m in range(len(multicol)):
                if m == 0:
                   s1 = ""; s2 = "|"; s3 = "|"
                elif m == len(multicol) - 1:
                   s1 = " & "; s2 = ""; s3 = "|"
                else:
                   s1 = " & "; s2 = ""; s3 = "||"
                s = s + "%s\\multicolumn{%d}{%sc%s}{%s}"%(s1,size[m],s2,s3,multicol[m])
            self.write(s)
            self.write("\\\\")
            self.write("\\hline")
        #else:
            #print "NO MULTICOL"

	first = False
        for mat in data:
            if first:
                self.write("\\hline")
            else:           
                first = True
                self.write("\\hline")
            for l in mat:
                if type(l[0]) == str:
                    ptype = "%s "
                elif isinstance(l[0],int):
                    ptype = "%d "
                else:
                    ptype = "%.2f "
                s = ptype%(l[0])
                #self.write("\\hline")

                for c in range(1,len(l)):
                    if type(l[c]) == str:
                        ptype = "& %s "
                    elif isinstance(l[c],int):
                        ptype = "& %d "
                    else:
                        ptype = "& %.2f "
                    #print printtype%(c)
                    s = s + ptype%(l[c])
                s = s + " \\\\"
                self.write(s)
                self.write("\\hline")

                      
        self.write("""
\\end{tabular}
\\caption{%s}
\\end{table}
"""%(caption))

    # image #
    def image(self,imagefile):
        self.write("\\begin{center}")
        self.write("\\includegraphics{%s}"%(imagefile))
        self.write("\\end{center}")

# close #
    def close(self):
        self.f.write(self.template[len(self.template) - 1])
        self.f.close()
        os.system("echo "" >> latex.output/tolatex.log")
        os.system("date >> latex.output/tolatex.log")
 
        # Open null file #
        FNULL = open(os.devnull, 'w')

        # Test if pdflatex is installed #
        rc = subprocess.call(['which', 'pdflatex'],stdout=FNULL, stderr=subprocess.STDOUT)
        if rc == 0:
	    if self.verbose:
		os.system("pdflatex -output-directory=. " + self.outputfile)
	    else:
		os.system("pdflatex -output-directory=. " + self.outputfile + " >> latex.output/tolatex.log")
	    #os.system("mv *.log *.aux *.tex latex.output"%(self.outputfile))
	    os.system("sleep 1")
	    os.system("mv *.log *.aux *.tex latex.output")
	    #os.system("mv latex.output/%s ."%(self.outputfile[:-4] + ".pdf"))
        else:
	    print("ERROR: It seems that pdflatex is not installed on your system.")
	    print("       Cannot export pdf files.")
	    print("       Please install pdflatex.")


######################################################################"

if __name__ == '__main__':
    
    tf = "|r|r|"
    mat = [[1,2],[3,4]]

    doc = latexdoc("test")
    doc.title("Test of the toprint class")
    doc.write("\\section{Hello}")
    doc.write("Bla bla bla $a_i$\\")
    doc.tableau([mat],tf,caption="test")
    doc.close()
