package indeldollo;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;


@Description("Tree likelihood for a cognate pruned tree")
public class PrunedTreeLikelihood extends TreeLikelihood {
	
	
	// Scaling parameters
	protected int X = 100;
	protected int m_nScale = 0;
	protected double m_fScale = 1.01;
	
	
	protected PrunedBeagleTreeLikelihood beaglePruned;
	
	protected CognatePrunedTree prunedTree;
	
	@Override
    public void initAndValidate() {
		
		
		// Ensure the tree is a cognate pruned tree
		if (!(treeInput.get() instanceof CognatePrunedTree)) {
			throw new IllegalArgumentException("Please ensure the tree is of type " + CognatePrunedTree.class.getCanonicalName());
		}
		this.prunedTree = (CognatePrunedTree) treeInput.get();
		
		
		if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
			String leaves = "?";
			if (treeInput.get() instanceof Tree) {
				leaves = String.join(", ", ((Tree) treeInput.get()).getTaxaNames());
			}
			throw new IllegalArgumentException(String.format(
					"The number of leaves in the tree (%d) does not match the number of sequences (%d). "
							+ "The tree has leaves [%s], while the data refers to taxa [%s].",
					treeInput.get().getLeafNodeCount(), dataInput.get().getTaxonCount(),
					leaves, String.join(", ", dataInput.get().getTaxaNames())));
		}
		
		beagle = null;
        beaglePruned = null;
        
     
        beaglePruned = new PrunedBeagleTreeLikelihood();
        try {
	        beaglePruned.initByName(
                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(), 
                    "useTipLikelihoods", m_useTipLikelihoods.get(), "scaling", scaling.get().toString(),
                    "rootFrequencies", rootFrequenciesInput.get());
	        if (beaglePruned.getBeagle() != null) {
	            //a Beagle instance was found, so we use it
	            return;
	        }
        } catch (Exception e) {
			// ignore
		}

        // No Beagle instance was found, so we use the good old java likelihood core
        beaglePruned = null;
        
        

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        //int stateCount = dataInput.get().getMaxStateCount();
        int stateCount = dataInput.get().getDataType().getStateCount();
        int patterns = dataInput.get().getPatternCount();
        likelihoodCore = createLikelihoodCore(stateCount);

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
		
	}
	
	
    @Override
	public double [] getPatternLogLikelihoods() {
		if (beaglePruned != null) {
			return beaglePruned.getPatternLogLikelihoods();
		}
		return patternLogLikelihoods.clone();
	} // getPatternLogLikelihoods

	
	
	 @Override
	 protected void initCore() {
	        final int nodeCount = prunedTree.getNodeCount();
	        
	        likelihoodCore.initialize(
	                nodeCount,
	                dataInput.get().getPatternCount(),
	                m_siteModel.getCategoryCount(),
	                true, m_useAmbiguities.get()
	        );

	        final int extNodeCount = nodeCount / 2 + 1;
	        final int intNodeCount = nodeCount / 2;

	        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
	            setPartials(prunedTree.getTrueRoot(), dataInput.get().getPatternCount());
	        } else {
	            setStates(prunedTree.getTrueRoot(), dataInput.get().getPatternCount());
	        }
	        hasDirt = Tree.IS_FILTHY;
	        for (int i = 0; i < intNodeCount; i++) {
	            likelihoodCore.createNodePartials(extNodeCount + i);
	        }
	        
	        if (scaling.get() == Scaling.always) {
	        	likelihoodCore.setUseScaling(1.01);
	        }
	        
	    }

	
	
	

    @Override
    public double calculateLogP() {
    	
    	
    	if (beaglePruned != null) {
            logP = beaglePruned.calculateLogP();
            return logP;
        }
    	
        // Unlike the standard tree likelihood, traversal begins at the birth node and not the root
    	this.prunedTree = (CognatePrunedTree) treeInput.get();
    	
    	//final Node root = prunedTree.getBirthRoot();
    	final Node root = prunedTree.getCognateMRCA();
    	
        ((PrunedTreeLikelihoodCore)this.likelihoodCore).setRootNr(root.getNr());

        try {
        	
        	// Set all branch lengths to -1 if the substitution matrix etc changed, including those above the mrca node
        	if (hasDirt != Tree.IS_CLEAN) {
        		Arrays.fill(m_branchLengths, -1);
        	}
        	
        	
        	// Negative branch length check
    		if (root.getHeight() > prunedTree.getBirthTime()) {
    			logP = Double.NEGATIVE_INFINITY;
    			return logP;
    		}
        	
        	
    		if (traverse(root, true) != Tree.IS_CLEAN) {
        		calcLogP();
        	}
        		
        }
        catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {

        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(root, true);
            calcLogP();
            return logP;
        }
        return logP;
    }
    
    


    /**
     * Traverse from the true root
     * Nodes outside of the cognate mrca do not contribute to the likelihood, however they still need to be checked for dirtiness for future iterations
     * @param node
     * @param isBirthRoot
     * @param followsBirthEvent
     * @return
     */
    protected int traverse(final Node node, boolean isBirthRoot) {
    	
    	

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        
        if (!node.isRoot() || isBirthRoot) {
        
	        final double parentHeight = isBirthRoot ? prunedTree.getBirthTime() : node.getParent().getHeight();
	        final double branchRate = branchRateModel.getRateForBranch(node);
	        final double branchTime = (parentHeight - node.getHeight()) * branchRate;
	
	        // First update the transition probability matrix(ices) for this branch
	        if (update  != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex]) {
	            
	        	update |= Tree.IS_DIRTY;
	        	
	            	
	            	m_branchLengths[nodeIndex] = branchTime;
	                likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
	            	
		            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
		                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
		                substitutionModel.getTransitionProbabilities(node, parentHeight, node.getHeight(), jointBranchRate, probabilities);
		                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
		            }
	            
	        }
        
        }
        
        

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {
        	
            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1, false);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2, false);
            
            
            // If either child node was updated then update this node too
            if (update != Tree.IS_CLEAN || update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }
                
                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                
                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                }

                if (isBirthRoot) {
                	
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    if (this.getConstantPattern() != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : this.getConstantPattern()) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                    double[] rootFrequencies = substitutionModel.getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }
                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                }

            }else {
 
        		
            }
        }
        return update;
    } 
    
    
    public void forceRecalculation() {
    	if (beaglePruned != null) {
            beaglePruned.forceRecalculation();
        }
    	hasDirt = Tree.IS_DIRTY;
    }
    
    
    @Override
	protected boolean requiresRecalculation() {
    	
    	if (beaglePruned != null) {
            return beaglePruned.requiresRecalculation();
        }
    	
    	boolean requiresRecalc = super.requiresRecalculation();
    	if (prunedTree.somethingIsDirty()) { 
			requiresRecalc = true;
		}
  
    	if (InputUtil.isDirty(prunedTree.birthInput)) {
    		requiresRecalc = true;
    	}

    	return requiresRecalc;
    }
    

	
	@Override
    protected LikelihoodCore createLikelihoodCore(int stateCount) {
		
		if (stateCount == 20) {
			return new PrunedTreeLikelihoodCore20();
		} else {
			return new PrunedTreeLikelihoodCore(stateCount);
		}
		
    }
	
	
	@Override
	public void store() {
		if (beaglePruned != null) {
            beaglePruned.store();
            //super.store();
            return;
        }
		super.store();
	}

	@Override
	public void restore() {
		if (beaglePruned != null) {
            beaglePruned.restore();
            //super.restore();
            return;
        }
		super.restore();
	}
}










