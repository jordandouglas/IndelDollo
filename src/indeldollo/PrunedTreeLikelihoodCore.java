package indeldollo;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.likelihood.BeerLikelihoodCore;



@Description("Tree likelihood core for a cognate pruned tree")
public class PrunedTreeLikelihoodCore extends BeerLikelihoodCore {
	
	
	boolean hasOrigin = true;
	int rootNr;
	int storedRootNr;
	
	double[] patternLogLikelihoodsOrigin;
	double[] storedPatternLogLikelihoodsOrigin;
	
	public PrunedTreeLikelihoodCore(int nrOfStates) {
		super(nrOfStates);
	}

	@Override
	public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities) {
		super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);
		
		this.patternLogLikelihoodsOrigin = new double[patternCount*nrOfStates];
		this.storedPatternLogLikelihoodsOrigin = new double[patternCount*nrOfStates];
		this.rootNr = -1;
		
	}
	
	
	public void setHasOrigin(boolean hasOrigin) {
		this.hasOrigin = hasOrigin;
	}
	
	
	/**
	 * Specify the pruned tree root number
	 * @param rootNr
	 */
	public void setRootNr(int rootNr) {
		this.rootNr = rootNr;
	}
	
	
	/**
	 * Calculate log-likelihoods, including the origin branch
	 */
	@Override
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        
		
        
        double[] matrix = matrices[currentMatrixIndex[rootNr]][rootNr];
        
        
        
        if (hasOrigin) {
        
	        double partialOriginBranch;
	
	        int u = 0;
	        int v = 0;
	
	        
	        // Calculate all partials at origin
	        for (int l = 0; l < 1; l++) { //nrOfMatrices; l++) {
	
	            for (int k = 0; k < nrOfPatterns; k++) {
	
	                int w = l * matrixSize;
	
	                for (int i = 0; i < nrOfStates; i++) {
	
	                	partialOriginBranch = 0.0;
	                    for (int j = 0; j < nrOfStates; j++) {
	                    	partialOriginBranch += matrix[w] * partials[v + j];
	                        w++;
	                    }
	
	                    this.patternLogLikelihoodsOrigin[u] = partialOriginBranch;
	                    u++;
	                }
	                v += nrOfStates;
	            }
	        }
	        
        
        }
        
        
        // Combine partials at origin
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            double sum = 0.0;
            for (int i = 0; i < nrOfStates; i++) {

            	//sum += frequencies[i] * partials[v];
                sum += frequencies[i] * (hasOrigin ? this.patternLogLikelihoodsOrigin[v] : 1);
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
    }
	
	

	
	
	
	 @Override
	 public void restore() {
		 double[] tmp = patternLogLikelihoodsOrigin;
		 patternLogLikelihoodsOrigin = storedPatternLogLikelihoodsOrigin;
		 storedPatternLogLikelihoodsOrigin = tmp;
		 
		 
		 //int tmp2 = rootNr;
		 rootNr = storedRootNr;
		 //storedRootNr = tmp2;
		 
		 super.restore();
	 }
	 
	 @Override
	 public void store() {
		 System.arraycopy(patternLogLikelihoodsOrigin, 0, storedPatternLogLikelihoodsOrigin, 0, patternLogLikelihoodsOrigin.length);
		 storedRootNr = rootNr;
		 super.store();
	 }


	 
	
	


}
