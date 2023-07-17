package indeldollo;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.core.Input.Validate;

@Description("Logs the genetic distance of each leaf to the root")
public class DistanceFromRootLogger extends BEASTObject implements Loggable {
	
	
	final public Input<Tree> treeInput = new Input<>("tree", "tree whose taxa distances to log", Validate.REQUIRED);
	final public Input<BranchRateModel> clockModelInput = new Input<>("clock", "clock model for the tree", Validate.OPTIONAL);

	
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void init(PrintStream out) {
		for (int leafNr = 0; leafNr < treeInput.get().getLeafNodeCount(); leafNr ++) {
			Node leaf = treeInput.get().getNode(leafNr);
			out.print("dist." + leaf.getID() + "\t");
		}
	}
	
	@Override
    public void log(long sample, PrintStream out) {
		for (int leafNr = 0; leafNr < treeInput.get().getLeafNodeCount(); leafNr ++) {
			Node leaf = treeInput.get().getNode(leafNr);
			out.print(getDistanceToRoot(leaf) + "\t");
		}
	}
	
	private double getDistanceToRoot(Node node) {
		
		if (node.isRoot()) {
			return 0;
		}
		
		double branchRate = clockModelInput.get() == null ? 1 : clockModelInput.get().getRateForBranch(node);
		double branchLength = node.getLength();
		double d = branchRate * branchLength;
		return d + getDistanceToRoot(node.getParent());
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
