package indeldollo.tools;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import indeldollo.util.PrunedTreeSet;

@Description("Reconstruct source tree from set of pruned trees")
public class SourceTreeFromPrunedTrees extends Runnable {
	final public Input<List<TreeFile>> treesInput = new Input<>("tree", "list of pruned tree files", new ArrayList<>(), Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified");
	
	final static double ROOT_HEIGHT = 1e10;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		
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
		Tree tree = treeSets.get(0).next();
		PrintStream out = System.out;
		if (outputInput.get() != null) {
			out = new PrintStream(outputInput.get());
		}
		tree.init(out);
		treeSets.get(0).reset();
		
        // recover taxon set
		List<Taxon> taxa = new ArrayList<>();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			String taxon = tree.getNode(i).getID();
			taxa.add(new Taxon(taxon));
		}
		TaxonSet taxonset = new TaxonSet(taxa);
		
		// merge pruned trees into source tree
		long k = 0;
        while (treeSets.get(0).hasNext()) {
        	// start with a caterpillar tree
    		Tree sourceTree = new Tree();
    		sourceTree.initByName("taxonset", taxonset);
    		for (int i = sourceTree.getLeafNodeCount(); i < sourceTree.getNodeCount(); i++) {
    			sourceTree.getNode(i).setHeight(ROOT_HEIGHT+i);
    		}

    		// add information from pruned trees one by one
    		for (int i = 0; i < treeSets.size(); i++) {
        		Tree prunedTree = treeSets.get(i).next();
    			merge(sourceTree, prunedTree);
    		}
    		
    		// log the result
    		sourceTree.log(k, out);
    		k++;
        }
        
		tree.close(out);
        
	}
	
	private void merge(Tree sourceTree, Tree tree) {		
		// sort nodes of tree to be merged in height
		Node [] nodes = tree.getNodesAsArray();
		int n = tree.getLeafNodeCount();
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
			Node node = nodes[i];
			if (node.getChildCount() > 0) {
				merge(sourceTree, node);
			}
		}
	}

	private void merge(Tree targetTree, Node node) {
		List<Node> prunedLeafs = node.getAllLeafNodes();
		List<Node> leafs = new ArrayList<>();
		Node [] nodes = targetTree.getNodesAsArray();
		for (Node n : prunedLeafs) {
			String id = n.getID();
			for (int i = 0; i < nodes.length; i++) {
				if (id.equals(nodes[i].getID())) {
					leafs.add(nodes[i]);
					break;
				}
			}
		}
		
		// find nodes to merge in target tree
		Node mrca = getMRCA(targetTree, leafs);
		if (mrca.getHeight() >= ROOT_HEIGHT) {
			mrca.setHeight(node.getHeight());
			while (mrca.getLeft().getHeight() >= ROOT_HEIGHT) {
				Node left = mrca.getLeft();
				if (!left.isLeaf()) {
					Node c = left.getLeft();
					if (c.isLeaf() && left.getRight().isLeaf() && !leafs.contains(c)) {
						c = left.getRight();
					}
					Node p = mrca.getParent();
	
					if (mrca.isRoot()) { // then p == null
						// p.removeChild(mrca);
						mrca.removeChild(left);
						left.removeChild(c);
						
						// p.addChild(left);
						mrca.addChild(c);
						left.addChild(mrca);
		
						c.setParent(mrca);
						mrca.setParent(left);
						
						left.setParent(null);
						targetTree.setRoot(left);
					} else {
						p.removeChild(mrca);
						mrca.removeChild(left);
						left.removeChild(c);
						
						p.addChild(left);
						mrca.addChild(c);
						left.addChild(mrca);
		
						c.setParent(mrca);
						left.setParent(p);
						mrca.setParent(left);
					}
				}
			}
		}
	}
	
    boolean [] nodesTraversed;
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
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	                nseen += 1;
	            }
	        } else if( h2 < h1 ) {
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
                    n = n1 = n.getParent();
                } else {
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
        nseen = 0;
        Node cur = leafs.get(0);

        for (int k = 1; k < leafs.size(); ++k) {
            cur = getCommonAncestor(cur, leafs.get(k));
        }
		return cur;
	}

	public static void main(String[] args) throws Exception {
		new Application(new SourceTreeFromPrunedTrees(), "", args);
	}
}
