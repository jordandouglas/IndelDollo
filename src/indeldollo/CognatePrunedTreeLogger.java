package indeldollo;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.TreeWithMetaDataLogger;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;


@Description("Logs a cognate pruned tree, but prunes above the cognate birth event, and below sampled death events")
public class CognatePrunedTreeLogger extends TreeWithMetaDataLogger {
	
	
	enum LetterCase{
		upper, lower, none
	}
	
	
	final public Input<Alignment> dataInput = new Input<>("data", "optional metadata to log at leaves", Validate.OPTIONAL);
	final public Input<String> removeInput = new Input<>("remove", "regex characters to remove from sequence (eg. gap characters)", Validate.OPTIONAL);
	final public Input<LetterCase> caseInput = new Input<>("case", "convert all sequences (if data != null) to upper case, lower case, or leave as they are (default)", 
															LetterCase.none, LetterCase.values());
	
    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;
    private DecimalFormat df;
    private boolean sortTree;
    
    
    final int maxNSampleAttempts = 1000;

    
    BranchRateModel clockModel;
    double mutationRate;
    
    
	 @Override
	 public void initAndValidate() {
		int dp = decimalPlacesInput.get();
        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }

        if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
        	someMetaDataNeedsLogging = false;
        }else {
        	someMetaDataNeedsLogging = true;
        }
    	
    	// Without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }

        // default is to sort the tree
        sortTree = false; //sortTreeInput.get();
        
        if ( !(treeInput.get() instanceof CognatePrunedTree)) {
        	throw new IllegalArgumentException("Please ensure that tree is a CognatePrunedTree");
        }
        
	        
	 }
	 
	 
	 
	 @Override
	 public void init(PrintStream out) {
		 CognatePrunedTree tree = (CognatePrunedTree) treeInput.get();
		 tree.getFullTree().init(out);
	 }
	
	
	@Override
    public void log(long sample, PrintStream out) {
		
        // make sure we get the current version of the inputs
		CognatePrunedTree tree = (CognatePrunedTree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + sample + " = ");

        if (sortTree) {
            tree.getRoot().sort();
        }
        
        
        // Start at the first birth event (not the root)
        Node node = tree.getBirthRoot();
        double birthTime = tree.getBirthTime();
        
        
        // Sample deaths and prune tree?
        List<Integer> nodesToPrune = new ArrayList<>();
        List<Double> timesToPrune = new ArrayList<>();


        
        out.print(toNewick(node, birthTime, metadata, branchRateModel, nodesToPrune, timesToPrune, true));
        out.print(";");
    }
	
	
	
	

    public String toNewick(Node node, double birthTime, List<Function> metadataList, BranchRateModel.Base branchRateModel, List<Integer> nodesToPrune, List<Double> timesToPrune, boolean isRoot) {
        
    	
    	// Prune thus subtree?
    	boolean prune = nodesToPrune.contains(node.getNr());
    	double pruneTime = -1;
    	if (prune) pruneTime = timesToPrune.get(nodesToPrune.indexOf(node.getNr()));
    	
    	StringBuffer buf = new StringBuffer();
        if (!prune && node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), birthTime, metadataList, branchRateModel, nodesToPrune, timesToPrune, false));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), birthTime, metadataList, branchRateModel, nodesToPrune, timesToPrune, false));
            }
            buf.append(")");
        } else {
        	
        	// If and this is not a leaf, then pick a leaf arbitrarily
        	int nodeNr = node.getNr();
        	if (prune && !node.isLeaf()) {
        		nodeNr = node.getAllLeafNodes().get(0).getNr();
        	}
            buf.append(nodeNr + 1);
        }
		StringBuffer buf2 = new StringBuffer();
		
		
		// Append sequence?
		String seq = "";
		if (node.isLeaf() && dataInput.get() != null && dataInput.get().getTaxonIndex(node.getID()) != -1) {
			
			// Get sequence
			seq = dataInput.get().getSequenceAsString(node.getID());
			
			// Filter out some characters
			if (removeInput.get() != null && !removeInput.get().isEmpty()) {
				seq = seq.replaceAll(removeInput.get(), "");
			}
			
			// Case conversion
			if (caseInput.get() == LetterCase.upper) {
				seq = seq.toUpperCase();
			}else if (caseInput.get() == LetterCase.lower) {
				seq = seq.toLowerCase();
			}
			
		}
		
		
		if (someMetaDataNeedsLogging || !seq.isEmpty()) {
			buf2.append("[&");
			if (metadataList.size() > 0) {
				for (Function metadata : metadataList) {
					if (metadata instanceof Parameter<?>) {
						Parameter<?> p = (Parameter<?>) metadata;
						int dim = p.getMinorDimension1();
						if (p.getMinorDimension2() > node.getNr()) {
							buf2.append(((BEASTObject) metadata).getID());
							buf2.append('=');
							if (dim > 1) {
								buf2.append('{');
								for (int i = 0; i < dim; i++) {
									if (metadata instanceof RealParameter) {
										RealParameter rp = (RealParameter) metadata;
										appendDouble(buf2, rp.getMatrixValue(node.getNr(), i));
									} else {
										buf2.append(p.getMatrixValue(node.getNr(), i));
									}
									if (i < dim - 1) {
										buf2.append(',');
									}
								}
								buf2.append('}');
							} else {
								if (metadata instanceof RealParameter) {
									RealParameter rp = (RealParameter) metadata;
									appendDouble(buf2, rp.getArrayValue(node.getNr()));
								} else {
									buf2.append(metadata.getArrayValue(node.getNr()));
								}
							}
						} else {
						
						}
					} else {
						if (metadata.getDimension() > node.getNr()) {
							buf2.append(((BEASTObject) metadata).getID());
							buf2.append('=');
							buf2.append(metadata.getArrayValue(node.getNr()));
						}
					}
					if (buf2.length() > 2 && metadataList.indexOf(metadata) < metadataList.size() - 1) {
						buf2.append(",");
					}
				}
				if (buf2.length() > 2 && branchRateModel != null) {
					buf2.append(",");
				}
			}
			

			if (!seq.isEmpty()) {
				buf2.append("seq='" + seq + "'");
				if (branchRateModel != null) buf2.append(",");
			}
			
			
			if (branchRateModel != null) {
				buf2.append("rate=");
				appendDouble(buf2, branchRateModel.getRateForBranch(node));
			}
			buf2.append(']');
		}
		if (buf2.length() > 3) {
			buf.append(buf2.toString());
		}
        
        //if (!isRoot) {
        	buf.append(":");
        	//double length = node.getLength();
        	double length = isRoot ? birthTime - node.getHeight() : node.getLength();
        	if (prune) length = node.getParent().getHeight() - pruneTime;
            if (substitutions) {
                appendDouble(buf, length * branchRateModel.getRateForBranch(node));
            } else {
                appendDouble(buf, length);
            }
        //}
       
        
        return buf.toString();
    }

    
    public void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }
	

}
