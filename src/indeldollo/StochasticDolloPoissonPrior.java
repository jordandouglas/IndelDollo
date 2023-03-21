package indeldollo;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.State;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;


@Description("Normalisation factor for stochastic dollo model, which is a function of the mean number of 1's expected in leaves according to a Poisson ")
public class StochasticDolloPoissonPrior extends TreeDistribution implements Loggable {
	
	public Input<RealParameter> lambdaInput = new Input<>("lambda", "cognate birth rate", Validate.REQUIRED);
	public Input<RealParameter> muInput = new Input<>("mu", "cognate death rate", Validate.REQUIRED);
	public Input<RealParameter> mutationRateInput = new Input<>("mutationRate", "to multiply by both rates", Validate.REQUIRED);
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel", "A model describing the rates on the branches of the beast.tree.");

    
    
    boolean needsUpdate = true;
    
    // Partials for u calculation
    boolean [] cladeDirty;
    
    
    protected double[] m_branchDeltas;
    protected double[] storedBranchDeltas;
    
    double[] uPartials;
    double[] stored_uPartials;
	
	
	@Override
	public void initAndValidate() {
		
		this.needsUpdate = true;
		if (treeInput.get() instanceof CognatePrunedTree) {
			throw new IllegalArgumentException("Please provide the full tree not the cognate pruned tree");
		}
		
		
		int nodeCount = this.treeInput.get().getNodeCount();
		this.uPartials = new double[nodeCount];
		this.stored_uPartials = new double[nodeCount];
		this.cladeDirty = new boolean[nodeCount];
		this.m_branchDeltas = new double[nodeCount];
		this.storedBranchDeltas = new double[nodeCount];
	}

	
	@Override
	public double calculateLogP() {
		
		
		if (this.needsUpdate) {
			
			
			
			double lambda = lambdaInput.get().getValue();
			double mu = muInput.get().getValue();
			Tree tree = (Tree) treeInput.get();
			BranchRateModel branchRateModel = branchRateModelInput.get();
			double mutationRate = mutationRateInput.get().getArrayValue();
			boolean dirtyRates = InputUtil.isDirty(muInput) || InputUtil.isDirty(mutationRateInput) || branchRateModelInput.get().isDirtyCalculation();
			this.calculateCladeDirtiness(tree.getRoot(), branchRateModel, mutationRate, mu, dirtyRates);
			
			double sum = 0;
			for (Node node_i : tree.getNodesAsArray()) {
				
				double x = 1; // If node is the root, then tj is infinity, hence x=1
				if (!node_i.isRoot()) {
					x = 1 - m_branchDeltas[node_i.getNr()];
				}
				x = x * (1-this.u(node_i));
				sum = sum + x; 
			}
			logP = -lambda/mu * sum;
			this.needsUpdate = false;
		}
		return logP;
		
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
		for (Node child : node.getChildren()) {
			double delta = Math.exp(-mu * mutationRate * child.getLength() * br.getRateForBranch(child));
			if (m_branchDeltas[child.getNr()] != delta) {
				dirty = true;
				m_branchDeltas[child.getNr()] = delta;
			}
		}
		
		
		this.cladeDirty[node.getNr()] = dirty;
		return this.cladeDirty[node.getNr()];
		
		
		
	}
	
	

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
	
	@Override
    public void store() {
        super.store();
        System.arraycopy(m_branchDeltas, 0, storedBranchDeltas, 0, m_branchDeltas.length);
        System.arraycopy(uPartials, 0, stored_uPartials, 0, uPartials.length);
        
    }

	 
	 
	
	@Override
    public void restore() {
		this.needsUpdate = true;
    	super.restore();
    	
    	double[] tmp = m_branchDeltas;
        m_branchDeltas = storedBranchDeltas;
        storedBranchDeltas = tmp;
        
        tmp = uPartials;
        uPartials = stored_uPartials;
        stored_uPartials = tmp;
        
	}
	
	@Override
    protected boolean requiresRecalculation() {
		
		//Log.warning(this.getID() + " requiresRecalculation");
		
		//this.needsUpdate = true;
        if (InputUtil.isDirty(mutationRateInput) || InputUtil.isDirty(muInput) || InputUtil.isDirty(lambdaInput) || treeInput.get().somethingIsDirty()) {
        	this.needsUpdate = true;
            return true;
        }

        if (branchRateModelInput.get().isDirtyCalculation()) {
        	this.needsUpdate = true;
            return true;
        }
        
        return this.needsUpdate;
	}
	
	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		return args;
	}


	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(mutationRateInput.get().getID());
		conds.add(treeInput.get().getID());
		conds.add(lambdaInput.get().getID());
		conds.add(muInput.get().getID());
		return conds;
	}
	
	@Override
	public void init(PrintStream out) {
		out.print(this.getID() + ".prob\t" + this.getID() + ".mean\t");
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		out.print(logP + "\t" + -logP + "\t");
	}
	

}






