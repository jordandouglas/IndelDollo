package simba;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;

@Description("Logs a cognate pruned tree's birth time")
public class CognateBirthLogger extends BEASTObject implements Loggable {
	
	
	final public Input<CognatePrunedTree> treeInput = new Input<>("tree", "cognate pruned tree whose birth time is logged", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void init(PrintStream out) {
		out.print(this.getID() + "\t");
	}
	
	@Override
    public void log(long sample, PrintStream out) {
		CognatePrunedTree tree = treeInput.get();
		double height = tree.getBirthRoot().getHeight();
		out.print(height + "\t");
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
