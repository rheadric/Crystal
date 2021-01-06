"""
    module Crystal
    export Atom, Cubic, FCC111, Deposit, PlotCryst
    export IsSupported, FunnelDown, Neighborhood

Crystal structures for simple cubic (001) and face-centered cubic (111) surfaces.
Includes functions to build a crystal by depositing atom-by-atom and to produce a plot.  

# Examples
```jldoctest
julia> using Crystal
[ Info: Precompiling Crystal [top-level]

julia> a = Atom(1,1,3,20,[])
Atom(1, 1, 3, 20, Process[])
 
julia> x = FCC111(3, 3)
FCC111(3, 3, Atom[Atom(1, 1, 1, 1, Process[#undef])], [1 1 1; 1 1 1; 1 1 1]

[0 0 0; 0 0 0; 0 0 0]

[0 0 0; 0 0 0; 0 0 0], 2)

julia> y = Cubic(4, 4)
Cubic(4, 4, Atom[Atom(1, 1, 1, 1, Process[#undef])], [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]

[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]

[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]

[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], 2)

julia> y.size
4

julia> y.layers
4

julia> y.AtomRegistry
1-element Array{Atom,1}:
 Atom(1, 1, 1, 1, Process[#undef])

 julia> y.world
 4×4×4 Array{Int64,3}:
 [:, :, 1] =
  1  1  1  1
  1  1  1  1
  1  1  1  1
  1  1  1  1
 
 [:, :, 2] =
  0  0  0  0
  0  0  0  0
  0  0  0  0
  0  0  0  0
 
 [:, :, 3] =
  0  0  0  0
  0  0  0  0
  0  0  0  0
  0  0  0  0
 
 [:, :, 4] =
  0  0  0  0
  0  0  0  0
  0  0  0  0
  0  0  0  0

 julia> y.nextatomID
2

julia> Deposit(y)
Atom(3, 3, 2, 2, Process[#undef])

julia> y.AtomRegistry
2-element Array{Atom,1}:
 Atom(1, 1, 1, 1, Process[#undef])
 Atom(3, 3, 2, 2, Process[#undef])
```
"""
module Crystal
export Atom, Cubic, FCC111, Deposit, PlotCryst, Process
export IsSupported, FunnelDown, Neighborhood
# 
# Abstract type AbstractCrystal with associated subtypes and functions.
# my goal is to produce a module that will support several different
# crystal types in a Kinetic Monte Carlo simulation.
# Following https://techytok.com/from-zero-to-julia/
# RLH 12/13/20
#

#
# Type AbstractCrystal is the base class for different crystal symmetries.
# The simplest concrete types are Cubic and Hexagonal.
# Other types are possible, such as FCC and CcCl.
# I plan to implement just two: Cubic and FCC(111).
# V0.1 RLH 12/15/20
# V0.2 RLH 12/19/20 Plan: Simple plotting for Cubic & FCC111
# V0.3 RLH 12/20/20 Plan: Implement IsSupported and FunnelDown
# V0.4 RLH 12/22/20 Implement Neighborhood() and AtomRegistry
# V0.5 RLH 12/23/20 Turn this into a module that I can load.
# V0.6 RLH 12/26/20 Atom now holds an array of type 'Process'.
#                   'Process' will be used by functions in the Processes module.
#                   A new verison of PlotCryst() takes an Atom as a second
#                   argument and makes a plot of just the neighborhood around
#                   that atom.

using PyPlot

"""
Holds process infotmation: AtomID, movevector, initialNN, finalNN, rate.

## Summary

mutable struct Process <: Any

## Fields

* AtomID     :: Int64
* movevector :: Array{Int64,1}
* initialNN  :: Int64
* finalNN    :: Int64
* rate       :: Float64
"""
mutable struct Process
    AtomID::Int64
    movevector::Array{Int64,1}
    initialNN::Int64
    finalNN::Int64
    rate:: Float64
end

# RLH 12/3/21 The code below helps splice!() to work on Arrays of type Process.
# I developed it without much understanding based on information found at:
# https://julialang.org/blog/2018/07/iterators-in-julia-0.7/
# The basic problem appears to be that splice!() is trying to iterate on the
# elements of the array rather than treating the as scalars.  I found several
# recommendations for how to turn off this behaviors but so far none have worked,
# aside from this method.
function Base.iterate(iter::Process, state=(1, 0)) 
    element, count = state
    if count > 1
        return nothing
    end
    return (iter[element], (iter[element], count + 1))
end
#
Base.length(y::Process) = 1
Base.getindex(y::Process, i::Int64) = y
Base.getindex(y::Process, z::Process) = z
#------------------------------------------------------------------------#

"""
A Struct holding the position and AtomID of an Atom.
Atoms are created by the Deposit() function or by calling the Atom() constructor.

# Summary
* mutable struct Atom <: Any
Fields
* xpos   :: Int64
* ypos   :: Int64
* zpos   :: Int64
* AtomID :: Int64
* processlist::Array{Process,1}

# Examples
```jldoctest
julia> julia> x = Cubic(3, 2)
Cubic(3, 2, Atom[Atom(1, 1, 1, 1, Process[#undef])], [1 1 1; 1 1 1; 1 1 1]

[0 0 0; 0 0 0; 0 0 0], 2)

julia> Deposit(x)
Atom(1, 1, 2, 2, Process[#undef])

Deposit(x)
Atom(1, 2, 2, 3, Process[#undef])

julia> x.world
3×3×2 Array{Int64,3}:
[:, :, 1] =
 1  1  1
 1  1  1
 1  1  1

[:, :, 2] =
 2  3  0
 0  0  0
 0  0  0

julia> x.AtomRegistry
3-element Array{Atom,1}:
 Atom(1, 1, 1, 1, Process[#undef])
 Atom(1, 1, 2, 2, Process[#undef])
 Atom(1, 2, 2, 3, Process[#undef])
```
"""
mutable struct Atom
    xpos::Int64
    ypos::Int64
    zpos::Int64
    AtomID::Int64
    processlist::Array{Process,1}
end

"""

### The base type for crystals.

Summary
* abstract type AbstractCrystal <: Any
    Subtypes 
* Cubic
* FCC111

"""
abstract type AbstractCrystal 
end

"""
A concrete subtype of AbstractCrystal.
Creates an (001) face of a simple cubic type crystal. 
Required arguments are size and layers.

# Summary
* mutable struct Cubic <: AbstractCrystal

Fields
* size       :: Int64
* layers     :: Int64
* AtomRegistry :: Array{Any,1}
* world      :: Array{Int64,3}
* nextatomID :: Int64

Supertype Hierarchy
* Cubic <: AbstractCrystal <: Any

# Notes
* Layer 1 is filled with immovable atoms by default and they all have AtomID 1. 
* The first deposited atom will land in layer 2 and it will have AtomID 2.

# Example
```jldoctest
julia> using Crystal
[ Info: Precompiling Crystal [top-level]

julia> x = Cubic( 10, 4)
Cubic(10, 4, Atom[Atom(1, 1, 1, 1, Process[#undef])], [1 1 … 1 1; 1 1 … 1 1; … ; 1 1 … 1 1; 1 1 … 1 1]

[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]

[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]

[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0], 2)

julia> for i=1:20 Deposit(x) end

julia> myfig = PlotCryst(x)
fig = Figure(PyObject <Figure size 1000x1000 with 1 Axes>)
```
"""
mutable struct Cubic <: AbstractCrystal 
    size::Int64
    layers::Int64
    AtomRegistry::Array{Atom,1}
    world::Array{Int64,3}
    nextatomID ::Int64
    # Constructor for Cubic
    function Cubic(size,layers)
        world = zeros(Int64,size,size,layers)
        nextatomID = 2
        # fill the first layer
        for i in 1:size
            for j in 1:size
                world[i,j,1] = 1
            end
        end
        AtomRegistry = [Atom(1,1,1,1,Array{Process}(undef, 1))]
        # The struct is created by new().
        new(size, layers, AtomRegistry, world, nextatomID)
    end
end

"""
A concrete subtype of AbstractCrystal.
It's just the same as Cubic, cuz I'm lazy.
Creates a face-centered cubic (111) type crystal. Required arguments are size and layers.

# Summary
* mutable struct FCC111 <: AbstractCrystal

Fields
* size       :: Int64
* layers     :: Int64
* AtomRegistry :: Array{Any,1}
* world      :: Array{Int64,3}
* nextatomID :: Int64

Supertype Hierarchy
* FCC111 <: Crystal <: Any

# Notes
* Layer 1 is filled with immovable atoms by default and they all have AtomID 1. 
* The first deposited atom will land in layer 2 and it will have AtomID 2.

# Example
```jldoctest
julia> using Crystal
[ Info: Precompiling Crystal [top-level]

julia> x = FCC111( 10, 4)
FCC111(10, 4, Atom[Atom(1, 1, 1, 1, Process[#undef])], [1 1 … 1 1; 1 1 … 1 1; … ; 1 1 … 1 1; 1 1 … 1 1]

[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]

[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]

[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0], 2)

julia> for i=1:20 Deposit(x) end

julia> myfig = PlotCryst(x)
fig = Figure(PyObject <Figure size 1000x1000 with 1 Axes>)
```
"""
mutable struct FCC111 <: AbstractCrystal
    size::Int64
    layers::Int64
    AtomRegistry::Array{Atom,1}
    world::Array{Int64,3}
    nextatomID ::Int64
    function FCC111(size, layers)
        x = Cubic(size, layers)
        new(x.size, x.layers, x.AtomRegistry, x.world, x.nextatomID)
    end
end


"""
Deposit(x::AbstractCrystal)

Deposit an atom onto crystal x at a random position.
This function accepts any crystal type as it's argument.
See the documentation for Cubic and FCC111 types for examples.
"""
function Deposit(x::AbstractCrystal)
    # choose a random position
    xpos, ypos = rand(1:x.size,2)
    # drop it onto the surface
    height =  x.layers
    while x.world[xpos, ypos, height] == 0
        height -= 1
    end 
    if height == x.layers
        error("Cannot exceed the maximum crystal height.")
    else
        atomID = x.nextatomID 
        x.world[xpos, ypos, height+1] = atomID
        procl = fill(Process(0, [0,0,0], 0, 0, 0.0),1);popfirst!(procl)
        newAtom = Atom(xpos, ypos, height+1, atomID, procl)
        push!(x.AtomRegistry, newAtom)
        x.nextatomID  += 1
    end
    return newAtom
end

"""
function PlotCryst(cubcryst::Cubic; labels = true, axis = true, topmax = 0)

Use scatter() to produce a simple plot of a Cubic crystal.
 """
function PlotCryst(cubcryst::Cubic; labels = true, axis = true, topmax = 0)
    #
    fig = figure("cubic_crystal",figsize=(10,10))
    ax = PyPlot.axes()
    #
    size = cubcryst.size
    areas = 200000.0/size^2
    #
    if topmax <= 0
        topmax = cubcryst.layers
    else
        topmax = min(topmax, cubcryst.layers)
    end
    #
    for k = 1:topmax
        x = []
        y = []
        for i = 1:size
            for j = 1:size
                if cubcryst.world[i,j,k] != 0
                    push!(x, i)
                    push!(y, j)
                end
            end
        end
        scatter(x,y,s=areas,alpha=1.0)
    end
    for i = 1:size
        for j = 1:size
            for k = topmax:-1:1
                if cubcryst.world[i,j,k] != 0
                    if labels 
                        ann = string(cubcryst.world[i,j,k])
                        annotate(ann,xy=[i-0.05;j-0.05])
                    end
                    break
                end
            end
        end
    end
    #
    PyPlot.title("Cubic Crystal")
    xlim((0.5,size+0.5))
    ylim((0.5,size+0.5))
    xlabel("X")
    ylabel("Y")
    if !axis
        grid("off")
        ax[:set_axis_off]()
    else
        grid("on")
    end
    @show fig
    return fig
end

"""
function PlotCryst(cubcryst::Cubic, myAtom::Atom; labels = true, axis = true)

Use scatter() to produce a simple plot of a Cubic crystal in the neighborhood around a given 'Atom'.
"""
function PlotCryst(cubcryst::Cubic, myAtom::Atom; labels = true, axis = true)
    #
    fig = figure("cubic_crystal_neighborhood",figsize=(10,10))
    ax = PyPlot.axes()
    #
    neighbors = convert(Array{Int64,3} ,  Neighborhood(cubcryst, myAtom))
    xpos = myAtom.xpos
    ypos = myAtom.ypos 
    zpos = myAtom.zpos
    print("myatom:", xpos," ", ypos," ",zpos,"\n")
    areas = 200000.0/3^2
    #
    for k = 1:3
        x = []
        y = []
        for i = 1:3
            for j = 1:3
                if neighbors[i,j,k] != 0
                    neighborAtom = cubcryst.AtomRegistry[neighbors[i,j,k]]
                    push!(x, xpos+i-2)
                    push!(y, ypos+j-2)
                end
            end
        end
        scatter(x,y,s=areas,alpha=1.0)
    end
    for i = 1:3
        for j = 1:3
            for k = 3:-1:1
                if neighbors[i,j,k] != 0
                    print(labels," ")
                    if labels 
                        neighborAtom = cubcryst.AtomRegistry[neighbors[i,j,k]]
                        print(i," ",j," ", k," ", neighborAtom.AtomID,"\n")
                        ann = string(neighborAtom.AtomID)
                        annotate(ann,xy=[xpos+i-2-0.05;ypos+j-2-0.05])
                    end
                    break
                end
            end
        end
    end
    #
    mystr = string("myatom:", xpos," ", ypos," ",zpos)
    PyPlot.title(mystr)
    xlim((xpos-1.5,xpos+1.5))
    ylim((ypos-1.5,ypos+1.5))
    xlabel("X")
    ylabel("Y")
    if !axis
        grid("off")
        ax[:set_axis_off]()
    else
        grid("on")
    end
    @show fig
    return fig
end

"""
IsSupported(x::FCC111,newAtom::Atom)

Check atoms below to see if the atom site is supported.
For FCC111 there may be up to three supporing atoms.

# Example
```jldoctest
julia> x=FCC111(3,2)
FCC111(3, 2, Atom[Atom(1, 1, 1, 1, Process[#undef])], [1 1 1; 1 1 1; 1 1 1]

[0 0 0; 0 0 0; 0 0 0], 2)

julia> a = Deposit(x)
Atom(2, 2, 2, 2, Process[#undef])

julia> IsSupported(x, a)
3
```
"""
function IsSupported(x::FCC111,newAtom::Atom)
    i = newAtom.xpos
    j = newAtom.ypos
    k = newAtom.zpos
    im = i<=1 ? x.size : i-1
    jm = j<=1 ? x.size : j-1
    ip = i>=x.size ? 1 : i+1
    jp = j>=x.size ? 1 : j+1
    supports = 0
    if mod(k,3) == 2
        # B layer i,j is supported by A layer
        # i,j i+1,j i+1,j+1
        #print(i," ",ip," ",j," ",jp," ",k,"\n")
        supports += x.world[i, j, k-1] > 0 ? 1 : 0
        supports += x.world[ip,j, k-1] > 0 ? 1 : 0
        supports += x.world[ip,jp,k-1] > 0 ? 1 : 0
    elseif mod(k,3) == 0
        # C layer i,j is supported by B layer
        # i,j  i-1,j i,j+1
        #print(i," ",im," ",j," ",jm," ",k,"\n")
        supports += x.world[i, j, k-1] > 0 ? 1 : 0
        supports += x.world[im,j, k-1] > 0 ? 1 : 0
        supports += x.world[i ,jp,k-1] > 0 ? 1 : 0
        #print(supports,"\n")
    else
        # A layer i,j is supported by C layer
        # i,j i-1,j-1 i,j-1
        supports += x.world[i, j, k-1] > 0 ? 1 : 0
        supports += x.world[im,jm,k-1] > 0 ? 1 : 0
        supports += x.world[i, jm,k-1] > 0 ? 1 : 0
    end
    return supports
end

"""
function FunnelDown(x::FCC111,newAtom::Atom)

Relax to a lower layer if atom is not supported.
For FCC111 there may be up to three supporing atoms. 
If there are fewer than three, the atom will funnel down to the next layer.

# Example
```julia
x = FCC111( 5, 5 )
for i=1:25
    myAtom = Deposit(x) 
    while IsSupported(x, myAtom) != 3
        FunnelDown(x, myAtom)
    end
end
```
"""
function FunnelDown(x::FCC111,newAtom::Atom)
    i = newAtom.xpos
    j = newAtom.ypos
    k = newAtom.zpos
    # If k=1 there is nothing to be done.
    if k==1 return -1 end
    # function wrap() enforces boundary conditions.
    wrap(i) = i>x.size ? 1 : (i<1 ? x.size : i)
    # Each layer type A,B,C is supported by the layer beneath.
    if mod(k,3) == 2
        iv = [0, 1, 1]
        jv = [0, 0, 1]
    elseif mod(k,3) == 0
        iv = [0, -1, 0]
        jv = [0,  0, 1]
    else
        iv = [0, -1,  0]
        jv = [0, -1, -1]
    end
    # Check IsSupported() to avoid possible infinite loop.
    if IsSupported(x,newAtom) == 3 return -3 end
    # Pick a site to funnel down to.
    while x.world[i, j, k] != 0
        ri = rand(1:3)
        inew = wrap(i+iv[ri])
        jnew = wrap(j+jv[ri])
        knew = k-1
        #print(inew," ",jnew," ",knew,"\n")
        if x.world[inew, jnew, knew] == 0
            x.world[inew, jnew, knew] = newAtom.AtomID
            x.world[i,j,k] = 0
            newAtom.xpos = inew
            newAtom.ypos = jnew
            newAtom.zpos = knew
        end
        #print(inew," ",jnew," ",knew,"\n")
    end
    # Recursive call if atom can funnel down again.
    #if IsSupported(x,newAtom) != 3
    #    FunnelDown(x, newatom)
    #end
end

"""
Return a chunk of world[] around a given atom.
"""
function Neighborhood(x::FCC111, newAtom::Atom)
    i = newAtom.xpos
    j = newAtom.ypos
    k = newAtom.zpos
    # function wrap() enforces boundary conditions.
    wrap(i) = i>x.size ? 1 : (i<1 ? x.size : i)
    # Need room for 6 nearest neighbors in the same layer.
    # Neigbors[:,2:3] will be for the same layer.
    Neighbors = zeros(Int, 3,4)
    if k>1 
        # For layer below, use the vectors from FunnelDown()
        # Each layer type A,B,C is supported by the layer beneath.
        if mod(k,3) == 2  # A<-B (B above A)
            iv = [0, 1, 1]
            jv = [0, 0, 1]
        elseif mod(k,3) == 0 #B<-C (C above B)
            iv = [0, -1, 0]
            jv = [0,  0, 1]
         else                #C<-A (A above C)
            iv = [0, -1,  0]
            jv = [0, -1, -1]
        end
        for int in 1:3
            ii = wrap(i+iv[int])
            jj = wrap(j+jv[int])
            Neighbors[int,1] = x.world[ii,jj,k-1]
        end
    end
    # There are 6 in-plane neighbors.
    iv = [1, 0, 1, -1, 0, 1]
    jv = [0, 1, 1, 0, -1, -1]
    for int in 1:3
        ii = wrap(i+iv[int])
        jj = wrap(j+jv[int])
        Neighbors[int,2] = x.world[ii,jj,k]
        ii = wrap(i+iv[int+3])
        jj = wrap(j+jv[int+3])
        Neighbors[int,3] = x.world[ii,jj,k]
    end
    if k<x.layers
        # For layer above use negative of FunnelDown() vectors.
        # Each layer type A,B,C is supported by the layer beneath.
        if mod(k,3) == 1   # A->B (A beneath B)
            iv = -[0, 1, 1]
            jv = -[0, 0, 1]
        elseif mod(k,3) == 2 #B->C (B beneath C)
            iv = -[0, -1, 0]
            jv = -[0,  0, 1]
        else                 #C->A (C beneath A)
            iv = -[0, -1,  0]
            jv = -[0, -1, -1]
        end
        for int in 1:3
            ii = wrap(i+iv[int])
            jj = wrap(j+jv[int])
            Neighbors[int,4] = x.world[ii,jj,k+1]
        end
    end
    return Neighbors
end

"""
IsSupported(x::Cubic,newAtom::Atom)

Check atoms below to see if the atom site is supported.
For Cubic there may be up to five supporing atoms, although only one can be directly beneath.
"""
function IsSupported(x::Cubic,newAtom::Atom)
    i = newAtom.xpos
    j = newAtom.ypos
    k = newAtom.zpos
    supports = 0
    ip = i>=x.size ? 1 : i+1
    jp = j>=x.size ? 1 : j+1
    im = i<=1 ? x.size : i-1
    jm = j<=1 ? x.size : j-1
    #print(i," ",im," ",j," ",jm," ",k,"\n")
    supports += x.world[i, j, k-1] > 0 ? 1 : 0
    supports += x.world[im,j, k-1] > 0 ? 1 : 0
    supports += x.world[i ,jm,k-1] > 0 ? 1 : 0
    supports += x.world[ip,j, k-1] > 0 ? 1 : 0
    supports += x.world[i ,jp,k-1] > 0 ? 1 : 0
    return supports
end

"""
function FunnelDown(x::Cubic,newAtom::Atom)

Relax to a lower layer if atom is not supported.
For Cubic there may be up to 5 supporing atoms, but in this implemntation only 2 are required to prevent downward funnelling.
One of the supporting atoms must be directly underneath.
"""
function FunnelDown(x::Cubic,newAtom::Atom)
    i = newAtom.xpos
    j = newAtom.ypos
    k = newAtom.zpos
    # If k=1 there is nothing to be done.
    if k==1 return -1 end
    # If not supported directly underneath, drop down.
    if x.world[i, j, k-1] == 0 
        x.world[i, j, k-1] = newAtom.AtomID
        x.world[i,j,k] = 0
        newAtom.zpos = k-1
        return 0
    end
    # function wrap() enforces boundary conditions.
    wrap(i) = i>x.size ? 1 : (i<1 ? x.size : i)
    iv = [1, 0, -1, 0]
    jv = [0, 1, 0, -1]
    # Check IsSupported() to avoid possible infinite loop.
    # If there is only just the 1 atom beneath, then relax.
    if IsSupported(x,newAtom) >=2 return -2 end
    # Pick a site to funnel down to.
    while x.world[i, j, k] != 0
        ri = rand(1:4)
        inew = wrap(i+iv[ri])
        jnew = wrap(j+jv[ri])
        knew = k-1
        #print(inew," ",jnew," ",knew,"\n")
        if x.world[inew, jnew, knew] == 0
            x.world[inew, jnew, knew] = newAtom.AtomID
            x.world[i,j,k] = 0
            newAtom.xpos = inew
            newAtom.ypos = jnew
            newAtom.zpos = knew
        end
        #print(inew," ",jnew," ",knew,"\n")
    end
end

"""
Return a chunk of world[] around a given atom.
"""
function Neighborhood(x::Cubic, newAtom::Atom)
    i = newAtom.xpos
    j = newAtom.ypos
    k = newAtom.zpos
    # function wrap() enforces boundary conditions.
    wrap(i) = i>x.size ? 1 : (i<1 ? x.size : i)
    Neighbors = zeros(Int,3,3,3)
    ii = map(wrap,i-1:i+1)
    jj = map(wrap,j-1:j+1)
    if k>1 
        Neighbors[:,:,1] = x.world[ii,jj,k-1]
    end
    Neighbors[:,:,2] = x.world[ii,jj,k]
    if k<x.layers
        Neighbors[:,:,3] = x.world[ii,jj,k+1]
    end
    return Neighbors
end

"""
function PlotCryst(fcc111cryst::FCC111; labels = true, axis = true, topmax=0)

Use PyPlot scatter() to produce a simple plot of an FCC111 crystal.
"""
function PlotCryst(fcc111cryst::FCC111; labels = true, axis = true, topmax=0)
    #
    fig = figure("fcc111_crystal",figsize=(10,sqrt(3)*5.2))
    ax = PyPlot.axes()
    #
    size = fcc111cryst.size
    areas = 150000.0/size^2
    # 
    if topmax <= 0
        topmax = fcc111cryst.layers
    else
        topmax = min(topmax, x.layers)
    end
    #
    for k = 1:topmax
        x = []
        y = []
        for i = 1:size
            for j = 1:size
                if fcc111cryst.world[i,j,k] != 0
                    xpos = i-j/2.0+0.5
                    ypos = sqrt(3)*j/2
                    # B & C layers get shifted in FCC111
                    if mod(k,3) == 2
                        xpos += 0.5
                        ypos += 0.5/sqrt(3)
                    elseif mod(k,3) == 0
                        ypos += 1.0/sqrt(3) 
                    end
                    if xpos < 0.99  
                        xpos += size 
                    end
                    push!(x, xpos)
                    push!(y, ypos)
                    if labels 
                        ann = string(fcc111cryst.world[i,j,k])
                        annotate(ann,xy=[xpos-0.05;ypos-0.05])
                    end
                end
            end
        end
        scatter(x,y,s=areas,alpha=0.8)
    end
    #
    PyPlot.title("FCC(111) Crystal")
    xlim((0.5,size+1.0))
    ylim(0.41,sqrt(3)*size/2.0+0.43+0.58)
    xlabel("X")
    ylabel("Y")
    if !axis
        grid("off")
        ax[:set_axis_off]()
    else
        grid("on")
    end
    @show fig
    return fig
end
#
end # end of module Crystal
