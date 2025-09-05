library(optparse)

# parse the options
option_list = list(
  make_option(c('-c','--start.count'),type='integer',default=1,help="Start count in the simulation"),
  make_option(c('-i','--start.iteration'),type='integer',default=1,help="Start iteration in the simulation (only used for count==start.count")
)

opt_parser = OptionParser(option_list=option_list)
option = parse_args(opt_parser)

# accessing the inputs
start.count = option$start.count
start.iteration = option$start.iteration


# define additional parameters
n2 = 100
n3 = 1000

p = 200
q = 4
s = 3

link = 'binomial'


source('simulation.R')









