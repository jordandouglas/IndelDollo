package indeldollo.tools;

import java.io.PrintStream;
import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beast.base.parser.NexusParser;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import indeldollo.util.PrunedTreeSet;

@Description("Reconstruct source tree set from set of pruned tree sets")
public class SourceTreeFromPrunedTrees extends Runnable {
	final public Input<List<TreeFile>> treesInput = new Input<>("tree", "list of pruned tree files", new ArrayList<>(), Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified");
	final public Input<Double> rootheightInput = new Input<>("rootheight", "minimum height of the root -- used when pruned tree have no information about root height", 2.0);
	
	double ROOT_HEIGHT;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		ROOT_HEIGHT = rootheightInput.get();
		
        // Open trees
        List<PrunedTreeSet> treeSets = new ArrayList<>();
        try {
        	for (TreeFile treefile : treesInput.get()) {
        		PrunedTreeSet treeSet = new PrunedTreeSet(treefile.getPath(), 0);
        		treeSets.add(treeSet);
            	treeSet.reset();
        	}
        } catch (Exception e) {
        	e.printStackTrace();
        	throw new Exception("Error Parsing Input Tree: " + e.getMessage());
        }
        
        // initialise output log
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
		}

		// recover taxon set
		TaxonSet taxonset = getTaxonSet();
		Tree sourceTree = new Tree();
		sourceTree.initByName("taxonset", taxonset);
		sourceTree.init(out);
		
		
		// merge pruned trees into source tree
		long k = 0;
        while (treeSets.get(0).hasNext()) {
        	// start with a caterpillar tree
    		sourceTree = new Tree();
    		sourceTree.initByName("taxonset", taxonset);
    		for (int i = sourceTree.getLeafNodeCount(); i < sourceTree.getNodeCount(); i++) {
    			sourceTree.getNode(i).setHeight(ROOT_HEIGHT+i*ROOT_HEIGHT/1e10);
    		}

    		// add information from pruned trees one by one
    		for (int i = 0; i < treeSets.size(); i++) {
        		Tree prunedTree = treeSets.get(i).next();
    			merge(sourceTree, prunedTree);
    		}
    		
    		// log the result
    		out.println();
    		sourceTree.log(k, out);
    		k++;
    		if (k % 10 == 0) {
        		if (k % 100 == 0) {
        			Log.warning.print("|");
        		} else {
        			Log.warning.print(".");
        		}
    		}
        }
        
		out.println();
		sourceTree.close(out);
        
		Log.warning("\nDone");
	}
	
	private TaxonSet getTaxonSet() {
		NexusParser parser = new NexusParser();
		try {
			parser.parseFile(treesInput.get().get(0));
		} catch (Throwable e) {
			// ignore
		}
		List<Taxon> taxa = new ArrayList<>();
		for (String taxon : parser.taxa) {
			taxa.add(new Taxon(taxon));
		}
		TaxonSet taxonset = new TaxonSet(taxa);
		return taxonset;
	}

	private void merge(Tree sourceTree, Tree prunedTree) {		
		// sort nodes of tree to be merged in height
		Node [] nodes = prunedTree.getNodesAsArray();
		int n = prunedTree.getLeafNodeCount();
		Arrays.sort(nodes, new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				// leave leaf nodes in place
				if (o1.getNr() < n && o2.getNr() < n) {
					return 0;
				}
				if (o1.getNr() < n) {
					return -1;
				}
				if (o2.getNr() < n) {
					return 1;
				}
				if (o1.getHeight() > o2.getHeight()) {
					return 1;
				}
				if (o1.getHeight() < o2.getHeight()) {
					return -1;
				}
				return 0;
			}
		});
				
		// merge nodes
		for (int i = n; i < nodes.length; i++) {
			merge(sourceTree, nodes[i]);
		}
	}

	private void merge(Tree sourceTree, Node node) {
		//System.out.println(sourceTree.getRoot().toNewick());
		

		List<Node> prunedLeafs = node.getAllLeafNodes();
		List<Node> leafs = new ArrayList<>();
		Node [] nodes = sourceTree.getNodesAsArray();
		for (Node n : prunedLeafs) {
			String id = n.getID();
//System.out.println(id);
			for (int i = 0; i < nodes.length; i++) {
				if (id.equals(nodes[i].getID())) {
					leafs.add(nodes[i]);
					break;
				}
			}
		}
		
		if (leafs.size() < 2) {
			int h = 3;
			h++;
		}
		
		// find nodes to merge in target tree
		Node mrca = getMRCA(sourceTree, leafs);
		if (mrca.getHeight() >= ROOT_HEIGHT) {
			// we need to add new node to tree
			if (subRoot.size() != 2) {
				int h = 3;
				getMRCA(sourceTree, leafs);
				h++;
			}
			Node [] subRootArray = subRoot.toArray(new Node[] {});
			Node left = subRootArray[0];
			Node right = subRootArray[1];
			Node lp = left.getParent();
			Node rp = right.getParent();
			if (lp == rp) {
				// no need to change topology
				mrca.setHeight(node.getHeight());				
			} else {
				int k = sourceTree.getLeafNodeCount();
				while (nodes[k].getHeight() < ROOT_HEIGHT) {
					k++;
				}
				Node n = nodes[k];
				if (n == lp) {
					// swap lp.otherchild with right
					rp.removeChild(right);
					Node other = otherChild(n, left);
					n.removeChild(other);
					rp.addChild(other);
					n.addChild(right);
					other.setParent(rp);
					right.setParent(n);
				} else if (n == rp) {
					// swap rp.otherchild with left
					lp.removeChild(left);
					Node other = otherChild(n, right);
					n.removeChild(other);
					lp.addChild(other);
					n.addChild(left);
					other.setParent(lp);
					left.setParent(n);
				} else {
					Node nLeft = n.getLeft();
					Node nRight = n.getRight();
					// swap n.left with left
					lp.removeChild(left);
					n.removeChild(nLeft);
					lp.addChild(nLeft);
					n.addChild(left);
					nLeft.setParent(lp);
					left.setParent(n);
					// swap n.right with right
					rp.removeChild(right);
					n.removeChild(nRight);
					rp.addChild(nRight);
					n.addChild(right);
					nRight.setParent(rp);
					right.setParent(n);
				}
				n.setHeight(node.getHeight());
			}
		}
//		System.out.println(sourceTree.getRoot().toNewick());
	}
	
    private Node otherChild(Node n, Node child) {
		if (n.getLeft() == child) {
			return n.getRight();
		}
		return n.getLeft();
	}

	boolean [] nodesTraversed;
    Set<Node> subRoot;
    protected int nseen;

    protected Node getCommonAncestor(Node n1, Node n2) {
        // assert n1.getTree() == n2.getTree();
        if( ! nodesTraversed[n1.getNr()] ) {
            nodesTraversed[n1.getNr()] = true;
            nseen += 1;
        }
        if( ! nodesTraversed[n2.getNr()] ) {
            nodesTraversed[n2.getNr()] = true;
            nseen += 1;
        }
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	        	if (!n1.isRoot() && n1.getParent().getHeight() >= ROOT_HEIGHT) {
	        		if (n1.getHeight() < ROOT_HEIGHT) {
	        			subRoot.add(n1);
	        		}
	        	}
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	                nseen += 1;
	            }
	        } else if( h2 < h1 ) {
	        	if (!n2.isRoot() && n2.getParent().getHeight() >= ROOT_HEIGHT) {
	        		if (n2.getHeight() < ROOT_HEIGHT) {
	        			subRoot.add(n2);
	        		}
	        	}
	            n2 = n2.getParent();
	            if( ! nodesTraversed[n2.getNr()] ) {
	                nodesTraversed[n2.getNr()] = true;
	                nseen += 1;
	            }
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
		        	if (!n1.isRoot() && n1.getParent().getHeight() >= ROOT_HEIGHT) {
		        		if (n1.getHeight() < ROOT_HEIGHT) {
		        			subRoot.add(n1);
		        		}
		        	}
                    n = n1 = n.getParent();
                } else {
    	        	if (!n2.isRoot() && n2.getParent().getHeight() >= ROOT_HEIGHT) {
    	        		if (n2.getHeight() < ROOT_HEIGHT) {
    	        			subRoot.add(n2);
    	        		}
    	        	}
                    n = n2 = n.getParent();
                }
	            if( ! nodesTraversed[n.getNr()] ) {
	                nodesTraversed[n.getNr()] = true;
	                nseen += 1;
	            } 
	        }
        }
        return n1;
    }

	protected Node getMRCA(Tree tree, List<Node> leafs) {
        nodesTraversed = new boolean[tree.getRoot().getAllChildNodesAndSelf().size()];
        subRoot = new HashSet<>();
        nseen = 0;
        Node cur = leafs.get(0);

        for (int k = 1; k < leafs.size(); ++k) {
            cur = getCommonAncestor(cur, leafs.get(k));
        }
		return cur;
	}

	public static void main(String[] args) throws Exception {
		new Application(new SourceTreeFromPrunedTrees(), "Source Tree From Pruned Trees", args);
	}
}
