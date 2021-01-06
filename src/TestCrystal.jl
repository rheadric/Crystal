push!(LOAD_PATH, "/Users/randallheadrick/Documents/myjulia/Crystal")
##
LOAD_PATH
##
using Crystal, PyPlot
##
print("hello world!")
##
x = Cubic( 10, 4)
for i=1:20 Deposit(x) end
#close(myfig)
myfig = PlotCryst(x)
##
IsSupported(x, Atom(8,5,3,13,[]))
##
FunnelDown(x, Atom(8,5,3,13,[]))
close(myfig)
myfig = PlotCryst(x)
##
x.world
##
close(myfig)
myfig = PlotCryst(x)
##
x = FCC111( 5, 5 )
for i=1:25
    myAtom = Deposit(x) 
    while IsSupported(x, myAtom) != 3
        FunnelDown(x, myAtom)
    end
end
close(myfig)
myfig = PlotCryst(x; labels=true, axis=true, topmax=0)
##
x = Cubic( 5, 5 )
for i=1:20
    myAtom = Deposit(x) 
    if IsSupported(x, myAtom) != 5
        print(myAtom.xpos," ",myAtom.ypos," ",myAtom.zpos,"\n")
        #FunnelDown(x, myAtom)
        #print(myAtom.xpos," ",myAtom.ypos," ",myAtom.zpos,"\n")
    end
end
close(myfig)
myfig = PlotCryst(x; labels=true, axis=true, topmax=0)
##
x = Cubic(5,5)
for i=1:20 Deposit(x) end
myfig = PlotCryst(x)
##
close(myfig)
##
myfig = PlotCryst(x, x.AtomRegistry[21])
##
PlotCryst(x, Atom(3,4,1,99,[]))
##
