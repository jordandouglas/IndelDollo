package simba;

import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("A tree with a variable root, based on a variable origin")
public class CognatePrunedTree extends Tree {
	
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree to prune", Input.Validate.REQUIRED);
	final public Input<RealParameter> birthInput = new Input<>("birth", "the time between the cognate mrca and the birth origin", Input.Validate.REQUIRED);
	final public Input<Alignment> dataInput = new Input<>("data", "a single site of cognate data", Validate.REQUIRED);
	
	final public Input<Boolean> samplingInput = new Input<>("sampling", "flag to indicate that this is for sampling only, so certain errors can be suppressed", false);
	
	
	Tree fullTree;
	boolean[] cognateAtLeaves; // Which taxa have 1's and which have 0's ?
	
	boolean[] cognateInClade;
	boolean[] storedCognateInClade;
	boolean[] deathInClade;
	boolean[] storedDeathInClade;
	
	
	
	double birthTime;
	double storedBirthTime;
	
	
	int birthNodeX; // Node where birth event occurred
	int mrcaNodeNr, storedMRCANodeNr; // Node of mrca
	
	boolean initialising = true; // Allow the likelihood to inisialise with the full node count
	
	// A dummy leaf with length 0 to enable an origin
	String dummyTaxon = "";
	final int dummyTaxonIndex = 0;
	
	
	int storedBirthRootNodeNr = -1;
	int birthRootNodeNr = -1;
	
	boolean allZeroes;
	private boolean needsUpdate = true;
	

	public CognatePrunedTree() {
		
	}
	
	
	@Override
    public void initAndValidate() {
		
		this.needsUpdate = true;
		//this.birthNode = -1;
		this.mrcaNodeNr = -1;
		this.fullTree = treeInput.get();
		this.initialising = true;
		
		Alignment data = dataInput.get();
		
		
		if (data.getSiteCount() != 1) {
			throw new IllegalArgumentException("Please ensure there is only 1 cognate in " + data.getID());
		}
		
		
		// Which taxa are 0's
		boolean allZeros = true;
		boolean allOnes = true;
		int numberOfOnes = 0;
		this.cognateAtLeaves = new boolean[data.getTaxonCount()];
		this.cognateInClade = new boolean[data.getTaxonCount()*2-1];
		this.storedCognateInClade = new boolean[data.getTaxonCount()*2-1];
		this.deathInClade = new boolean[data.getTaxonCount()*2-1];
		this.storedDeathInClade = new boolean[data.getTaxonCount()*2-1];		

		for (int j = 0; j < data.getTaxonCount(); j ++) {
			String taxon = this.getNode(j).getID();
			String sequence = data.getSequenceAsString(taxon);
			//List<Integer> symbols = data.getDataType().stringToEncoding(sequence);
			//int symbolInt = symbols.get(0);
			
			Log.warning(taxon + " " + sequence + " for taxon " + (j+1) + " in " + data.getID());
			
			if (sequence.equals("0")) {
				this.cognateAtLeaves[j] = false;
				allOnes = false;
			}else {
				this.cognateAtLeaves[j] = true;
				allZeros = false;
				numberOfOnes ++;
			}
		}
		Log.warning(this.getID() + ": " + numberOfOnes + " of " + data.getTaxonCount() + " taxa have 1's");
		
		allZeroes = false;
		if (!samplingInput.get()) {
			if (allZeros) {
				allZeroes = true;
				Log.warning("Warning: there should be at least one '1' at the leaves");
			}
			
			if (allOnes) {
				Log.warning("Warning: could not find any '0's' at the leaves for " + this.getID());
			}
		}
		
		// Birth node of the cognate
		// Initialise it above the true root so that the starting tree has the full taxon set (to simplify the initialisation of other classes)
		//birthInput.get().setValue(this.fullTree.getRoot().getHeight() * 2);
		//this.birthNode = this.fullTree.getRoot().getNr();
		this.mrcaNodeNr = this.fullTree.getRoot().getNr();
		
		
		// Find a unique name for the dummy sequence
		List<String> taxonNames = Arrays.asList(this.fullTree.getTaxaNames());
		this.dummyTaxon = "dummy";
		while (taxonNames.contains(this.dummyTaxon)) {
			this.dummyTaxon += "0";
		}
		
		
	}
	
	

	
	 /**
	  * Get the node where the cognate was born above
	  * @param siteNum
	  * @return
	  */
	public int getBirthNodeNr() {
		
		this.updateBirthNodes();
		int nr = this.mrcaNodeNr;
		if (nr < 0) {
			Log.warning("warning: birthNode is -1");
			nr = this.fullTree.getRoot().getNr();
		}
		return nr;
	}
	
	@Override
	public void getMetaData(final Node node, final Double[] t, final String pattern) {
		this.fullTree = treeInput.get();
		this.fullTree.getMetaData(node, t, pattern);
	}
	
	

    
    @Override
    protected boolean requiresRecalculation() {
    	
    	if (isDirty()) {
    		this.needsUpdate = true;
    		return true;
    	}
    	return true;
    }
    
    @Override
    public void setEverythingDirty(final boolean isDirty) {
    	
    }
    
    
    @Override
	public boolean somethingIsDirty() {
		boolean d = isDirty();
		if (d) this.needsUpdate = true;
		return d;
	}
    
    private boolean isDirty() {
    	if (dataInput.get().isDirtyCalculation()) return true;
		if (InputUtil.isDirty(birthInput)) return true;
		
		
		/*
		// Recalculate the birth node
		int old = birthNode;
		double oldTime = birthTime;
		this.needsUpdate = true;
		this.updateBirthNodes();
		
		
		
		if (old != birthNode) dirty = true;
		if (oldTime != birthTime) dirty = true;
		birthNode = old;
		birthTime = oldTime;
		*/
		
		
		//boolean dirty = false;
		//if (this.fullTree.getNode(birthNode).isDirty() != Tree.IS_CLEAN) {
			//dirty = true;
		//}
		
		
		return this.fullTree.somethingIsDirty();
		//return cladeIsDirty(this.fullTree.getNode(birthNode));
    }
    
    private boolean cladeIsDirty(Node node) {
    	if (node.isDirty() != Tree.IS_CLEAN) return true;
    	if (node.isLeaf()) return false;
    	if (cladeIsDirty(node.getChild(0))) return true;
    	if (cladeIsDirty(node.getChild(1))) return true;
    	return false;
    }

    
    @Override
    public void setSomethingIsDirty(final boolean isDirty) {
    	//this.fullTree.setSomethingIsDirty(isDirty);
    	//super.setSomethingIsDirty(isDirty);
    	if (isDirty) this.needsUpdate = true;
    }
	 
    
    @Override
    public Tree copy() {
    	return this; //.fullTree.copy();
    }
    
    
    public Tree getFullTree() {
    	this.fullTree = treeInput.get();
    	return this.fullTree;
    }
    
    
    @Override
   	public String toString() {
    	return getRoot().toString();
    }
    
    
    /**
     * The root of the full tree
     * @return
     */
	public Node getTrueRoot() {
		this.fullTree = treeInput.get();
		return this.fullTree.getRoot();
	}
	
	

	
	/**
     * The root under the first birth event
     * @return
     */
	public Node getBirthRoot() {
		return this.fullTree.getNode(getBirthNodeNr());
	}
	
	/**
     * The mrca of all 1's
     * @return
     */
	public Node getCognateMRCA() {
		this.updateBirthNodes();
		return this.fullTree.getNode(this.mrcaNodeNr);
	}
	

    
	@Override
	public Node getRoot() {
		return this.getBirthRoot(); 
	}
	
	
	/**
	 * Get the node nr of this node in the full tree
	 * @return
	 */
//	public int getFullTreeNodeNr(int nr) {
//		if (nr == getBirthNodeNr()) return - 1; // New root
//		return nr;
//	}
	
	
	@Override
	public String [] getTaxaNames() {
		return this.fullTree.getTaxaNames();
	}

	 /**
	  * Find the node corresponding to the birth on this tree
	  */
	 private void updateBirthNodes() {
		 
		 this.fullTree = treeInput.get();
		 
		if (!this.needsUpdate) return;
		
		
		// Calculate clade cognates
		cladeHasCognate(this.fullTree.getRoot());
		cladeHasDeath(this.fullTree.getRoot());
		this.needsUpdate = false;
		
		 
		if (this.cognateAtLeaves == null) {
			//this.birthNode = -1;
			this.mrcaNodeNr = -1;
			return;
		}
		
		
		// Pathological case. Just use the root
		if (allZeroes) {
			//this.birthNode = 
			this.mrcaNodeNr = this.getFullTree().getRoot().getNr();
			return;
		}
		
		
		// Find the most recent ancestor such that all of the 1's are beneath it
		Node node = findBirthNode(this.fullTree.getRoot());
		
		if (node == null) {
			//this.birthNode = -1;
			this.mrcaNodeNr = -1;
		}else {
			
			this.mrcaNodeNr = node.getNr();
			
			// Go back to the first node after the birth event
			this.birthTime = birthInput.get().getValue() + node.getHeight();
			while (!node.isRoot() && node.getParent().getHeight() < birthTime) {
				node = node.getParent();
			}
			//this.birthNode = node.getNr();
		}
		
		
		
			 
		
	}
	
	
	/**
	 * Get MRCA of all 1's
	 * @param node
	 * @param siteNum
	 * @return
	 */
	public Node findBirthNode(Node node) {
		
		if (node.isLeaf()) {
			if (this.cognateAtLeaves[node.getNr()]) return node;
			return null;
		}
		
		
		node.getChild(0);
		node.getChild(1);
		
		boolean leftHasCognate = cladeHasCognate(node.getChild(0));
		boolean rightHasCognate = cladeHasCognate(node.getChild(1));
		
		
		// Found the node
		if (leftHasCognate && rightHasCognate) return node;
		
		if (leftHasCognate) return findBirthNode(node.getChild(0));
		if (rightHasCognate) return findBirthNode(node.getChild(1));
		
		// Pathological case: no cognate (all 0's)
		return null;
		
		
	}
	
	

	/**
	 * Does the clade have a 1?
	 **/
	public boolean cladeHasCognate(Node node) {
		if (!this.needsUpdate) {
			return cognateInClade[node.getNr()];
		}
		
		boolean val = false;
		if (node.isLeaf()) {
			val = this.cognateAtLeaves[node.getNr()];
		}else {
		
			if (cladeHasCognate(node.getChild(0))) {
				val = true;
			}
			
			if (cladeHasCognate(node.getChild(1))) {
				val = true;
			}
			
		}
		cognateInClade[node.getNr()] = val;
		return val;
	}
	
	
	/**
	 * Does the clade have a 0?
	 **/
	public boolean cladeHasDeath(Node node) {
		
		if (!this.needsUpdate) {
			return deathInClade[node.getNr()];
		}
		
		boolean val = false;
		if (node.isLeaf()) {
			val = !this.cognateAtLeaves[node.getNr()];
		}else {
			
			if(cladeHasDeath(node.getChild(0))) {
				val = true;
			} 
			
			if(cladeHasDeath(node.getChild(1))) {
				val = true;
			}
		}
		deathInClade[node.getNr()] = val;
		return val;
	}
	 

	
	@Override
	public void store() {
		
		//this.fullTree.setEverythingDirty(false);
		this.getFullTree().getRoot().makeAllDirty(Tree.IS_CLEAN);
		storedBirthRootNodeNr = birthRootNodeNr;
		storedMRCANodeNr = mrcaNodeNr;
		
		System.arraycopy(cognateInClade, 0, storedCognateInClade, 0, cognateInClade.length);
		System.arraycopy(deathInClade, 0, storedDeathInClade, 0, deathInClade.length);
		
		
		storedBirthTime = birthTime;
		
	}
	


	@Override
	public void restore() {
		
		hasStartedEditing = false;
		postCache = null;
		
		//this.fullTree.setEverythingDirty(false);
		this.getFullTree().getRoot().makeAllDirty(Tree.IS_CLEAN);
		
		//super.restore();
		birthRootNodeNr = storedBirthRootNodeNr;
		mrcaNodeNr = storedMRCANodeNr;
		//this.needsUpdate = true;
		
		boolean[] tmp = cognateInClade;
		cognateInClade = storedCognateInClade;
		storedCognateInClade = tmp;
		
		tmp = deathInClade;
		deathInClade = storedDeathInClade;
		storedDeathInClade = tmp;
		
		
		birthTime = storedBirthTime;
		

	}
	
	@Override
	public int getLeafNodeCount() {
		return this.fullTree.getLeafNodeCount(); 
	}
	
	@Override
	public int getInternalNodeCount() {
		return this.fullTree.getInternalNodeCount(); 
	}
	
	@Override
	public int getNodeCount() {
		return this.fullTree.getNodeCount(); 
	}
	
	

	
	@Override
	public Node getNode(int i) {
		return this.fullTree.getNode(i);
	}
	
	@Override
	public Node[] getNodesAsArray() {
		return this.fullTree.getNodesAsArray();
	}
		
	
	

	/**
	 * The time of the cognate birth
	 * @return
	 */
	public double getBirthTime() {
		this.updateBirthNodes();
		return this.birthTime;
	}
	
	
	/**
	 * The index of the dummy taxon
	 * @return
	 */
	public int getDummyNr() {
		return this.dummyTaxonIndex;
	}
	
	
	@Override
	public void fromXML(final org.w3c.dom.Node node) {
		// Do nothing
	}
	
	
	@Override
	public void assignFromFragile(final StateNode other) {
		// Do nothing
	}
	
	
	/**
	 * The dummy taxon name
	 * @return
	 */
	public String getDummyTaxon() {
		return this.dummyTaxon;
	}


	public void stopInit() {
		this.initialising = false;
	}


	public boolean allZeros() {
		return this.allZeroes;
	}
	


	

	
}




