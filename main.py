import EasySectionAnalysis as cs
import numpy as np
import matplotlib.pyplot as plt
import os

def general_plotter(
    plots, 
    title:str=None,
    xlabel:str=None,
    xlim:tuple=None,
    xinvert:bool=False,
    ylabel:str=None,
    ylim:tuple=None,
    yinvert:bool=False,
    grid:bool=False,
    legend=False, 
    fname:str=None, 
    show:bool=True
    ):
    """This function is a general plotting function.
    It sets up a plot with the given parameters.
    Compress plot setup down to one function. for more varies option use pyplot by itself

    Arguments:
        plots {list of tuples} -- [(xs,ys,label,(marker))]
    
    Keyword Arguments:
        title {str} -- title displayed above the plot (default: {None})
        xlabel {str} -- label of the x-axis (default: {None})
        xlim {tuple} -- the limit of the x-axis (default: {None})
        xinvert {bool} -- invert the x-axis (default: {False})
        ylabel {str} -- label of the y-axis(default: {None})
        ylim {tuple} -- the limit of the y-axis (default: {None})
        yinvert {bool} -- invert the y-axis (default: {False})
        grid {bool} -- turns on or off the grid (default: {False})
        legend {bool or int} -- turns on the legend if passed a value. An int maybe passed for locating the legend. This is the same as pyplot (default: {False})
        fname {str} -- filename to save an image of the plot. If passed it will try to save the file with that filename (default: {None})
        show {bool} -- Show the plot in it's own window. Set false if executing in seperate threads (default: {True})
    """
    for plot in plots:
        if len(plot) == 4:
            xs,ys,label,marker = plot
        else:
            xs,ys,label = plot
            marker="None"
        plt.plot(xs,ys,label=label,marker=marker)
    if title!=None:
        plt.title(title)
    if xlabel !=None:
        plt.xlabel(xlabel)
    if xlim!=None:
        plt.xlim(xlim)
    if xinvert:
        plt.gca().invert_xaxis()
    if ylabel!=None:
        plt.ylabel(ylabel)
    if ylim != None:
        plt.ylim(ylim)
    if yinvert:
        plt.gca().invert_yaxis()
    if grid:
        plt.grid()
    #setup legend
    if type(legend)==int:
        plt.legend(loc=legend)
    else:
        if legend:
            plt.legend(loc=0) 
    #save the figure with fname
    if fname!=None: 
        plt.savefig(fname)
    else:
        if not show:
            print("Why do you want to create a graph that you don't save or show.\nThis is utterly useless")
    if show: 
        plt.show()
    plt.close()

def latex_fig_nx2(graphs:list):
    imgstr = []
    for j in range(0,len(graphs)-(len(graphs)%6),6):
        imgstr.append("\\begin{figure}[t!] % \"[t!]\" placement specifier just for this example")
        for i in range(0,6,2):
            line = '''\\begin{{subfigure}}{{0.48\\textwidth}}
                \\includegraphics[width=\\linewidth]{{{0}}}
                \\caption{{cross-section of section {1}}} \\label{{fig:crossection{1}}}
                \\end{{subfigure}}\\hspace*{{\\fill}}
                \\begin{{subfigure}}{{0.48\\textwidth}}
                \\includegraphics[width=\\linewidth]{{{2}}}
                \\caption{{cross-section of section {3}}} \\label{{fig:crossection{3}}}
                \\end{{subfigure}}
                \\medskip'''.format(graphs[j+i],j+i+1,graphs[j+i+1],j+i+2)
            imgstr.append(line)
        imgstr.append("\\end{figure}")
    return "\n".join(imgstr)


def calculateWing(infile):
    sections = np.genfromtxt(infile)
    casename = infile.split(".")[0]
    saveFiles = True
    sectioncoef = []
    pltnames=[]
    sc_xs8=[]
    sc_ys8=[]
    sc_xs=[]
    sc_ys=[]
    for row in sections:
        # row:
        #   0   1      2   3    4     5     6   7    8           9       10
        # tskin,Ireq,Nsec,chord,t/c,Hfront, w, tspar,Astringer,N_str_top,N_str_bot,
        sec = cs.WinboxSection()
        sec.Nstringerstop = int(row[9]) # number of stringers on the upper panel of the wingbox
        sec.Nstringersbottom = int(row[10]) # number of stringers on the bottom panel of the wingbox
        sec.Astringertop = row[8] # [m^2] the cross sectional area of the top stringer
        sec.Astringerbottom = row[8] # [m^2] the cross sectional area of the bottom stringer
        sec.Asparcap = row[8] # [m^2] I think this one has become obselete
        sec.As_sparcaps = [row[8], # [m^2] the cross sectional area of top right sparcap
                            row[8], # [m^2] the cross sectional area of top left sparcap
                            row[8], # [m^2] the cross sectional area of bottom left
                            row[8]] # [m^2] the cross sectional area of bottom right
        sec.Hfront = row[5] # [m] height of the front spar
        sec.Hback = row[5]#*0.8 # [m] height of the back spar
        sec.tSpar = row[7] # [m] thickness of the spar
        sec.tskin = row[0] # [m] thickness of the skin
        sec.c = row[6] # legnth of the root chord

        sec.createGeometry()

        sectioncoef.append(sec.calculateParameters())
        pltnames.append("plot_section_{}.png".format(int(row[2])))

        c,I,SC,NA  = sectioncoef[-1]
        sc_x,sc_y = SC
        sc_xs8.append(sc_x/sec.c)
        sc_ys8.append((sec.Hfront/2+sc_y)/sec.Hfront)

        if saveFiles:
            sec.plotGeometry(fname="{0}/sections/{0}_plot_section_{1}.png".format(casename, int(row[2])))
            sec.saveBoomtoFile("{0}/sections/{0}_boom_section_{1}.csv".format(casename, int(row[2])))
            sec.saveParameterstoFile("{0}/sections/{0}_pars_section_{1}.csv".format(casename, int(row[2])))


    xs = np.linspace(0.51375, 41.1/2,num=20)
    general_plotter(
        [
            (xs,sections[:,1],"required stiffness","."),
            (xs,[x[1][0] for x in sectioncoef],"stiffness with adjusted N.A.","^")
    ],
    title="required stiffness and calculated stiffness along the half wing",
    xlabel="[m]",
    ylabel="Ixx [m^4]",
    legend=0,
    grid=True,
    fname="{0}/{0}_Ixx-plot.png".format(casename),
    show=False

    )

    # plt.scatter(sc_xs8,sc_ys8,label="shearcenters of each section",marker="o", c="b")
    # plt.scatter(sum(sc_xs8) / len(sc_xs8), sum(sc_ys8) / len(sc_ys8),label="average shearcenter",marker="x", c="r")
    # plt.title("position of shearcenter for back-spar = 0.8*front-spar")
    # plt.xlabel("x/c")
    # plt.ylabel("y/(t/c)")
    # plt.legend()
    # plt.grid()
    # plt.show()
    # plt.close()

def main():
    cases = [
        "Case1.csv",
        "Case2.csv",
        "Case3.csv",
        "Case4.csv",
        "Case5.csv"
    ]
    for case in cases:
        calculateWing(case)
    # print(latex_fig_nx2(["Images/crossections/{}".format(x) for x in pltnames]))
    # wingbox = section.WinboxSection()
    # wingbox.createGeometry()
    # wingbox.plotGeometry()
    
if __name__ == "__main__":
    main()