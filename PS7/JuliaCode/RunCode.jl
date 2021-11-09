
# theme(:juno)
# default(fontfamily="Computer Modern", framestyle=:box)
# Load the Code
    include("Code.jl");

#Run The Code
    #All exercises are done at once like this since they are supposed to all use
        #the same e draw
    #Note that the results are highly dependent upon the epsilon draw,
    #estpecially for exercise 6. Oftentimes standard errors cannot be
    #calculated. Setting UseRandomSeed to false gives a consistent draw where this
    #is not the case. Refresh Julia after you do this to reset the seed somewhere
    #more truly random
    StepsAThroughD(UseRandomSeed=false)
    StepsAThroughD(UseRandomSeed=true)
