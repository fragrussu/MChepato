### Generate vertices for the hepatocyte meshes from simple regular prism references



## import modules
import numpy as np
import sys

## input and output folders
reffolder = '../data/ref'
pertfolder = '../data/perturbed'

## seed for reproducibility
np.random.seed(19870210)

## cell sizes [um]
cellsize = np.concatenate(  ( np.linspace(5.0,50.0,31) , np.array([51.5,53.0,54.5,56.0,57.5,60.0]) )  )
nsides = np.array([4,5,6])

## number of random cell perturbations
npert = 5

## extent of the perturbation
fpert = 0.1  # perturbation will be fpert*cellsize


## Rescale cells and perturbe each size
print('')
print('+++++++++++++++++++++++++++++++++++++++++++')
print('        Generating mesh verices            ')
print('+++++++++++++++++++++++++++++++++++++++++++')
print('')

# loop through shapes (4, 5 and 6 sides)
for pp in range(0,nsides.size):

	print('    ... {} sides'.format(nsides[pp]))

	# load 1um reference
	refvert = np.loadtxt( '{}/sides{}_diam1um.vertices'.format(reffolder,nsides[pp]) )

	# loop through hepatocyte sizes
	for cc in range(0,cellsize.size):
	
		# get actual size in um
		mysize = float(cellsize[cc])
		print('             ... size {} um'.format(mysize))

		# rescale vertices
		myvert = mysize*refvert

		# store in the reference folder
		np.savetxt('{}/sides{}_diam{}um.vertices'.format(reffolder,nsides[pp],mysize), myvert, fmt='%.8f', delimiter=' ')

		# now perturb the regular vertices
		for nn in range(0,npert):

			print('                       ... perturbation {}'.format(nn))

			# perturb with normally distributed numbers
			myvertpert = myvert + (fpert*mysize)*np.random.randn(myvert.shape[0],myvert.shape[1])

			# store in the perturbation folder
			np.savetxt('{}/sides{}_diam{}um_fpert{}_npert{}.vertices'.format(pertfolder,nsides[pp],mysize,fpert,nn), myvertpert, fmt='%.8f', delimiter=' ')		



print('')
print('        Done')
print('')
print('')

