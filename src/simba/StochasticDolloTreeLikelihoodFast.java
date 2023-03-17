package simba;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;


/**
 * Based on
 * Nicholls, Geoff K., and Russell D. Gray. 
 * "Dated ancestral trees from binary trait data and their application to the diversification of languages." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70, no. 3 (2008): 545-566.
 */
@Description("Stochastic dollo tree likelihood with birth times explicitly sampled")
public class StochasticDolloTreeLikelihoodFast extends GenericTreeLikelihood {

	
    public Input<RealParameter> muInput = new Input<RealParameter>("mu", "cognate death rate", Validate.REQUIRED);
    public Input<RealParameter> lambdaInput = new Input<RealParameter>("lambda", "cognate birth rate", Validate.REQUIRED);
    public Input<RealParameter> mutationRateInput = new Input<RealParameter>("rate", "relative clock rate");
    public Input<RealParameter> cognateBirthInput = new Input<RealParameter>("birth", "cognate birth", Validate.OPTIONAL);
    public Input<CognatePrunedTree> cognateTreeInput = new Input<CognatePrunedTree>("cognateTree", "cognate pruned tree", Validate.REQUIRED);
    
    protected BranchRateModel.Base branchRateModel;
    
    boolean stochasticIsDirty;
    int zeroState = -1;
    int oneState = -1;
    
    CognatePrunedTree prunedTree;
    
    
    boolean debug = false; 
    
    
    // Using 2 doubles instead of double array of length 2 for efficiency
    double zeroPartials; 
    double onePartials;
    
    
    // Per node
    protected double[] nodePartialsOne;
    protected double[] storedNodePartialsOne;
    protected double[] nodePartialsZero;
    protected double[] storedNodePartialsZero;
    
    protected double[] m_branchDeltas;
    protected double[] storedBranchDeltas;
    
    // Partials for u calculation
    double[] uPartials;
    double[] stored_uPartials;
    
    boolean [] cladeDirty;
    
    
    public StochasticDolloTreeLikelihoodFast() {
    	treeInput.setRule(Validate.OPTIONAL);
    	siteModelInput.setRule(Validate.OPTIONAL);
    }
    
    
	@Override
	public void initAndValidate() {
		
		debug = false; // this.getID().equals("treeLikelihood.trait.ala");

		
		
		
		this.stochasticIsDirty = true;
		
		// Pruned tree
		this.prunedTree = cognateTreeInput.get();
		
		 // Sanity check: alignment should have same #taxa as tree
		if (dataInput.get().getTaxonCount() != this.prunedTree.getLeafNodeCount()) {
			String leaves = "?";
			if (this.prunedTree.getFullTree() instanceof Tree) {
				leaves = String.join(", ", this.prunedTree.getTaxaNames());
			}
			throw new IllegalArgumentException(String.format(
					"The number of leaves in the tree (%d) does not match the number of sequences (%d). "
							+ "The tree has leaves [%s], while the data refers to taxa [%s].",
					this.prunedTree.getFullTree().getLeafNodeCount(), dataInput.get().getTaxonCount(),
					leaves, String.join(", ", dataInput.get().getTaxaNames())));
		}
		
		
		if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
		
		// Ensure the datatype is binary
		DataType dt = dataInput.get().getDataType();
		if (dt.getStateCount() != 2) {
			throw new IllegalArgumentException("Please ensure the datatype is binary (there are " + dt.getStateCount() + " states but there should be 2)");
		}
		if (dt.getCharacter(0).equals("0")){
			this.zeroState = 0;
			this.oneState = 1;
		}else if (dt.getCharacter(1).equals("0")){
			this.zeroState = 1;
			this.oneState = 0;
		}else {
			throw new IllegalArgumentException("Please ensure one of the two symbols is a 0 (the current states are " + dt.getCharacter(0) + " and " + dt.getCharacter(1) + ")");
		}
		
		int nodeCount = this.prunedTree.getFullTree().getNodeCount();
		this.uPartials = new double[nodeCount];
		this.stored_uPartials = new double[nodeCount];
		
		this.nodePartialsOne = new double[nodeCount];
		this.storedNodePartialsOne = new double[nodeCount];	
		this.nodePartialsZero = new double[nodeCount];
		this.storedNodePartialsZero = new double[nodeCount];	
		
		this.cladeDirty = new boolean[nodeCount];
		this.m_branchDeltas = new double[nodeCount];
		this.storedBranchDeltas = new double[nodeCount];
		
	}

	
	
	@Override
	public double calculateLogP() {
		
		if (debug) Log.warning("calculateLogP");
		
		//if (!stochasticIsDirty) return logP;
		
		branchRateModel = branchRateModelInput.get();
		
		 
		final Node root = this.prunedTree.getRoot();
		 
		
		// Negative branch length check
		for (Node node : this.prunedTree.getNodesAsArray()) {
			double length = node.getLength();
			if (length < 0) {
				logP = Double.NEGATIVE_INFINITY;
				return logP;
			}
		}
		

		// Get parameters
		double lambda = lambdaInput.get().getValue();
		double mu = muInput.get().getValue();
		double birthTime = prunedTree.getBirthTime();
		double mutationRate = mutationRateInput.get() == null ? 1 : mutationRateInput.get().getArrayValue();
		
		// Check birth node is above the root 
		if (mutationRate <= 0 || lambda <= 0 || mu <= 0 || birthTime < root.getHeight()) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}
		
		
		boolean dirtyRates = prunedTree.somethingIsDirty() || InputUtil.isDirty(muInput) || InputUtil.isDirty(mutationRateInput);
		this.calculateCladeDirtiness(root, branchRateModel, mutationRate, mu, dirtyRates);
		
		
		// Felsenstein's from the root
		this.traverse(root, mu, mutationRate);
		logP = this.onePartials;
		
		
		
		// Mean birth rate (weighted by branch lengths)
        double totalWeightedRate = 0.0;
        double totalTreeLength = 0.0;
        for (Node node : root.getAllChildNodesAndSelf()) {
        	double rate = branchRateModel.getRateForBranch(node);
        	double length = node == root ? (birthTime - node.getHeight()) : node.getLength();
            totalWeightedRate += rate * length;
            totalTreeLength += length;
        }
        double meanBranchRate = totalWeightedRate/totalTreeLength;
        double meanBirthRate = lambda * mutationRate * meanBranchRate;
        
		//double meanBirthRate = lambda * mutationRate * branchRateModel.getRateForBranch(root);
		
		
		// Origin branch
		double branchTime = (birthTime - root.getHeight()) * mutationRate * branchRateModel.getRateForBranch(root) ;
		logP += -branchTime*mu;
		
		
		logP += Math.log(meanBirthRate); 
		
		// Ascertainment correction (i.e. remove probability of being all zeroes)
		logP += Math.log(1 - u(root)); 
	
				
		
		stochasticIsDirty = false;
		return logP;
		
	
	}
	


	
	
	protected void traverse(Node node, double deathRate, double clockRate) {
		
		
		boolean p = false; //this.getID().equals("treeLikelihood.cognates.2.5");
		
		if (!this.cladeDirty[node.getNr()]) {
			//this.onePartials = this.nodePartialsOne[node.getNr()];
			//this.zeroPartials = this.nodePartialsZero[node.getNr()];
			//return;
		}
		
		
		
		// Case 1: this is a leaf
		if (node.isLeaf()) {
			if (this.prunedTree.cladeHasCognate(node)) {
				this.onePartials = 0;
				this.zeroPartials = Double.NEGATIVE_INFINITY;
			}else {
				this.onePartials = Double.NEGATIVE_INFINITY;
				this.zeroPartials = 0;
			}
			if (p) Log.warning("case1: " + this.zeroPartials + ", " + this.onePartials);
			this.nodePartialsOne[node.getNr()] = this.onePartials;
			this.nodePartialsZero[node.getNr()] = this.zeroPartials;
			return;
		}
		
		
		
		
		
		
		// Case 2: there are no zeroes in the clade
		// Only need to look at the probability of not dying
		if (!this.prunedTree.cladeHasDeath(node)) {
			this.onePartials = traverseBirthsOnly(node, deathRate, clockRate);
			this.zeroPartials = Double.NEGATIVE_INFINITY;
			this.nodePartialsOne[node.getNr()] = this.onePartials;
			this.nodePartialsZero[node.getNr()] = this.zeroPartials;
			if (p) Log.warning("case2: " + this.onePartials);
			return;
		}
		
		
		
		// Case 3: there is a 0 and a 1 in the clade
		// Regular Felsenstein's likelihood algorithm, but hardcoded to 2 states, and incorporating knowledge of irreversible deaths
		double onePartialsThis = 0;
		double zeroPartialsThis = 0;
		double branchTime, logP11, logP10;;
		Node child;
		for (int i = 0; i < 2; i ++) {
			
			child = node.getChild(i);
			traverse(child, deathRate, clockRate);
			
			// 0 partials for this node. Redundant calculation commented out is shown here for completeness
			// double logP00 = 0;
			// double logP01 = Double.NEGATIVE_INFINITY;
			// zeroPartialsThis += this.logsumexp(logP00+this.zeroPartials, logP01+this.onePartials);
			zeroPartialsThis += this.zeroPartials;
			
			
			// 1 partials for this node.
			branchTime = child.getLength() * clockRate * branchRateModel.getRateForBranch(child);
			logP11 = -branchTime*deathRate;
			logP10 = Math.log(1 - Math.exp(-branchTime*deathRate));
			onePartialsThis += this.logsumexp(logP10+this.zeroPartials, logP11+this.onePartials);
			
			
			if (p) Log.warning(node.getNr() + " case3: " + child.getNr() + " - " +  + zeroPartialsThis + ", " + onePartialsThis);

			
		}
		
		this.onePartials = onePartialsThis;
		this.zeroPartials = zeroPartialsThis;
		
		this.nodePartialsOne[node.getNr()] = this.onePartials;
		this.nodePartialsZero[node.getNr()] = this.zeroPartials;

		
	}
	
	
	/**
	 * Returns log(exp(val1) + exp(val2))
	 */
	private double logsumexp(double val1, double val2) {
		if (val1 == Double.NEGATIVE_INFINITY) return val2;
		if (val2 == Double.NEGATIVE_INFINITY) return val1;
		return Math.log(1 + Math.exp(val2 - val1)) + val1;
	}
	
	
	/**
	 * Fast likelihood calculation conditional on there being no 0's in the clade
	 * Deaths are irreversible, and therefore we only need to look at branch lengths
	 */
	private double traverseBirthsOnly(Node node, double deathRate, double clockRate) {
		
		if (node.isLeaf()) return 0;
		double logPStayAlive = 0;
		for (int i = 0; i < 2; i ++) {
			Node child = node.getChild(i);
			double branchTime = child.getLength() * clockRate * branchRateModel.getRateForBranch(child);
			logPStayAlive += -branchTime*deathRate; // Log prob of not dying in time t with deathrate mu
			logPStayAlive += traverseBirthsOnly(child, deathRate, clockRate);
		}
		return logPStayAlive;
		
	}


	/**
	 * Recursively calculate the 'u0' term from equation (5)
	 * https://rss.onlinelibrary.wiley.com/doi/epdf/10.1111/j.1467-9868.2007.00648.x
	 */
	private double u(Node node) {
		
		
		if (node.isLeaf()) {
			return 0;
		}
		
		if (!this.cladeDirty[node.getNr()]) {
			return this.uPartials[node.getNr()];
		}
		
		double result = 1;
		for (int c = 0; c < 2; c++) {
			
			// Child calculation
			Node child = node.getChild(c);
			double delta = m_branchDeltas[child.getNr()];
			double r = 1 - delta + delta*u(child);
			result = result * r;
		}
		
		
		this.uPartials[node.getNr()] = result;
		this.cladeDirty[node.getNr()] = false;
		return result;
		
		
	}
	
	
	/**
	 * Flag all nodes as dirty if they or one of their descendants are dirty
	 * 
	 * @param node
	 * @return
	 */
	private boolean calculateCladeDirtiness(Node node, BranchRateModel br, double mutationRate, double mu, boolean dirtyRates) {

		
		// Leaves are always clean
		if (node.isLeaf()) {
			this.cladeDirty[node.getNr()] = false;
			return false;
		}
		
		boolean dirty = dirtyRates;
		
		// Is left clade dirty?
		if (calculateCladeDirtiness(node.getChild(0), br, mutationRate, mu, dirtyRates)){
			dirty = true;
		}
		
		// Is right clade dirty?
		if (calculateCladeDirtiness(node.getChild(1), br, mutationRate, mu, dirtyRates)){
			dirty = true;
		}
		
		// Is this node dirty?
		if (node.isDirty() != Tree.IS_CLEAN) {
			dirty = true;
		}
		
		// Branch rates different?
		for (int c = 0; c < 2; c++) {
			Node child = node.getChild(c);
			double delta = Math.exp(-mu * mutationRate * child.getLength() * br.getRateForBranch(child));
			if (m_branchDeltas[child.getNr()] != delta) {
				dirty = true;
				m_branchDeltas[child.getNr()] = delta;
			}
			
		}
		
		
		this.cladeDirty[node.getNr()] = dirty;
		return dirty;
		
		
		
	}
	
	
	
	
	
	
	@Override
    protected boolean requiresRecalculation() {
		
		if (debug) Log.warning("requiresRecalculation");

		
        if (InputUtil.isDirty(lambdaInput)) {
        	stochasticIsDirty = true;
            return true;
        }
        if (InputUtil.isDirty(muInput)) {
        	stochasticIsDirty = true;
            return true;
        }
        if (mutationRateInput.get() != null && InputUtil.isDirty(mutationRateInput)) {
        	stochasticIsDirty = true;
            return true;
        }
        
        if (InputUtil.isDirty(branchRateModelInput)) {
        	stochasticIsDirty = true;
            return true;
        }
        
        if (prunedTree.somethingIsDirty()) { 
        	stochasticIsDirty = true;
            return true;
        }
        
        return stochasticIsDirty;
        
    }
	
	
	@Override
	public List<String> getArguments() {
		return super.getArguments();
	}


	@Override
	public List<String> getConditions() {
		List<String> conds = super.getConditions();
		conds.add(lambdaInput.get().getID());
		conds.add(muInput.get().getID());
		conds.add(prunedTree.birthInput.get().getID());
		return conds;
	}
	
	
	
	@Override
    public void accept() {
		if (debug) Log.warning("accept");
		super.accept();
	}
	
	@Override
    public void store() {
		if (debug) Log.warning("store");
        super.store();
        System.arraycopy(m_branchDeltas, 0, storedBranchDeltas, 0, m_branchDeltas.length);
        System.arraycopy(uPartials, 0, stored_uPartials, 0, uPartials.length);
        System.arraycopy(nodePartialsOne, 0, storedNodePartialsOne, 0, nodePartialsOne.length);
        System.arraycopy(nodePartialsZero, 0, storedNodePartialsZero, 0, nodePartialsZero.length);
    }

    @Override
    public void restore() {
    	if (debug) Log.warning("restore");
    	stochasticIsDirty = true;
        super.restore();
        
    	double[] tmp = m_branchDeltas;
        m_branchDeltas = storedBranchDeltas;
        storedBranchDeltas = tmp;
        
        tmp = uPartials;
        uPartials = stored_uPartials;
        stored_uPartials = tmp;
        
        tmp = nodePartialsOne;
        nodePartialsOne = storedNodePartialsOne;
        storedNodePartialsOne = tmp;
        
        tmp = nodePartialsZero;
        nodePartialsZero = storedNodePartialsZero;
        storedNodePartialsZero = tmp;

        
    }



	
	
}












