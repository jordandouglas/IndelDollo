package simba;

import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Node;

@Description("Initialises a birth time in cognate birth model (above the mrca of the 1's)")
public class CognateBirthModelInitialiser extends BEASTObject implements StateNodeInitialiser {

	
	final public Input<RealParameter> birthTimeInput = new Input<>("birthTime", "birth time of cognate", Input.Validate.OPTIONAL);
	final public Input<CognatePrunedTree> treeInput = new Input<>("tree", "cognate pruned tree", Input.Validate.OPTIONAL);
	
	CognatePrunedTree tree;
	RealParameter birthTime;
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public void initStateNodes() {
		
		// Get the pruned tree and birth time
		tree = treeInput.get();
		birthTime = birthTimeInput.get();
		
		
		// Find the mrca of the 1's
		Node mrca = tree.findBirthNode(tree.getFullTree().getRoot());
		
		
		// Set the birth time to slightly above the mrca height
		double treeHeight = tree.getFullTree().getRoot().getHeight();
		double time = (mrca == null || mrca.isRoot()) ? treeHeight * 1.01 : mrca.getHeight() +  Randomizer.nextFloat() * mrca.getLength();
		birthTime.setValue(time);
		
		Log.warning("Setting " + birthTime.getID() + " to height " + time);
		
		
		
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(birthTime);
	}



}
