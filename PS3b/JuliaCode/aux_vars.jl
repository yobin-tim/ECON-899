# Identifiers
indentifiers = Dict(:market => :Year, :product => :id)

# Name of covariates to use for the regression
covariates_names = [:Intercept
                    :delta_iia
                    :price
                    :dpm
                    :hp2wt
                    :size
                    :turbo
                    :trans
                    :Year_1986
                    :Year_1987
                    :Year_1988
                    :Year_1989
                    :Year_1990
                    :Year_1991
                    :Year_1992
                    :Year_1993
                    :Year_1994
                    :Year_1995
                    :Year_1996
                    :Year_1997
                    :Year_1998
                    :Year_1999
                    :Year_2000
                    :Year_2001
                    :Year_2002
                    :Year_2003
                    :Year_2004
                    :Year_2005
                    :Year_2006
                    :Year_2007
                    :Year_2008
                    :Year_2009
                    :Year_2010
                    :Year_2011
                    :Year_2012
                    :Year_2013
                    :Year_2014
                    :Year_2015
                    :model_class_2
                    :model_class_3
                    :model_class_4
                    :model_class_5
                    :cyl_2
                    :cyl_4
                    :cyl_6
                    :cyl_8
                    :drive_2
                    :drive_3]


# Name of the instruments
instruments_names = [:i_import
                    :diffiv_local_0
                    :diffiv_local_1
                    :diffiv_local_2
                    :diffiv_local_3
                    :diffiv_ed_0
                    :blpiv_0
                    :blpiv_1
                    :blpiv_2
                    :blpiv_3
                    :blpiv_4
                    :blpiv_5
]

sim_data_names = [:Var1]

model_specs = Dict( :ids         => indentifiers,
                    :covariates  => covariates_names,
                    :instruments => instruments_names,
                    :sim_data    => sim_data_names)
                   