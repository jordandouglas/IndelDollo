package indeldollo.tools;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import beast.base.util.CollectionUtils;
import beast.base.util.DiscreteStatistics;
import beast.base.util.HeapSort;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.CladeSystem;
import beastfx.app.treeannotator.ContourMaker;
import beastfx.app.treeannotator.ContourPath;
import beastfx.app.treeannotator.ContourWithSynder;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import indeldollo.util.PrunedTreeSet;

@Description("Finds maximum clade credibility tree among the most commonly-occurring taxon set of a pruned tree posterior")
public class PrunedTreeAnnotator extends beast.base.inference.Runnable {

	
	public Input<String> treesInput = new Input<>("trees", "NEXUS file containing a tree set", Input.Validate.REQUIRED);
	public Input<String> outputInput = new Input<>("out", "output file", Input.Validate.REQUIRED);
	public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);

	final boolean forceIntegerToDiscrete = false;
	final double hpd2D = 0.8;
    int totalTrees = 0;
    int totalTreesUsed = 0;
	Set<String> attributeNames = new HashSet<>();
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() throws Exception {
		

		CladeSystem cladeSystem = new CladeSystem();
		
        totalTrees = 10000;
        totalTreesUsed = 0;
        int burninPercentage = burnInPercentageInput.get();

        
        // Open trees
        PrunedTreeSet treeSet;
        try {
        	treeSet = new PrunedTreeSet(treesInput.get(), burninPercentage);
        } catch (Exception e) {
        	e.printStackTrace();
        	throw new Exception("Error Parsing Input Tree: " + e.getMessage());
        }
        
        
        
        // Each tree may have a different taxonset. Use the most commonly occurring one
        List<Boolean> includeTree = filterTreeset(treeSet);
        
        
        // Count number of trees
        try {
        	treeSet.reset();
            cladeSystem.setProcessSA(false);
            int treeNum = 0;
        	while (treeSet.hasNext()) {
        		
        		Tree tree = treeSet.next();
        		
        		// Include tree?
        		if (includeTree.get(treeNum)) {
	        		setupAttributes(tree);
	                tree.getLeafNodeCount();
	            	cladeSystem.add(tree, false);
	            	totalTreesUsed++;
        		}
        		treeNum++;
               
            }
            totalTrees = totalTreesUsed * 100 / (100-Math.max(burninPercentage, 0));
        } catch (Exception e) {
        	throw new Exception(e.getMessage());
        }

        if (totalTrees < 1) {
        	throw new Exception("No trees");
        }
        if (totalTreesUsed <= 1) {
            if (burninPercentage > 0) {
            	throw new Exception("No trees to use: burnin too high");
            }
        }
        
        
        

        // Clade probabilities
        cladeSystem.calculateCladeCredibilities(totalTreesUsed);
        System.out.println("Total number of trees " + totalTrees + ", where " + totalTreesUsed + " are used.");
        System.out.println("Total unique clades: " + cladeSystem.getCladeMap().keySet().size());
        System.out.println();

        
        // Get MCC tree and sort it
        System.out.println("Finding maximum credibility tree...");
        Tree targetTree = summarizeTrees(treeSet, cladeSystem, includeTree).copy();
        targetTree.getRoot().sort();
        
        
        
        System.out.println("Collecting node information...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");
        
        
        // Collect node information
        int stepSize = Math.max(totalTreesUsed / 60, 1);
        int reported = 0;

        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        cladeSystem = new CladeSystem();
        cladeSystem.setProcessSA(false);
        cladeSystem.add(targetTree, true);
        int totalTreesUsedNew = 0;
        try {
            int counter = 0;
            treeSet.reset();
            while (treeSet.hasNext()) {
            	Tree tree = treeSet.next();
            	
            	if (!includeTree.get(counter)) {
            		counter++;
            		continue;
            	}
                
                cladeSystem.collectAttributes(tree, attributeNames);
                if (counter > 0 && counter % stepSize == 0 && reported < 61) {
    				while (1000 * reported < 61000 * (counter + 1)/ this.totalTreesUsed) {
    					System.out.print("*");
    	                reported++;
            	    }
    				System.out.flush();
                }
                totalTreesUsedNew++;
                counter++;
        	}
        	
            cladeSystem.removeClades(targetTree.getRoot(), true);
            this.totalTreesUsed = totalTreesUsedNew;
            cladeSystem.calculateCladeCredibilities(totalTreesUsedNew);
        } catch (Exception e) {
            Log.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        System.out.println();
        System.out.println();

        System.out.println("Annotating target tree...");

        
        // Annotate the tree
        try {
            annotateTree(cladeSystem, targetTree.getRoot(), null);
            setTreeHeightsByCA(treeSet, targetTree, includeTree);
        } catch (Exception e) {
        	e.printStackTrace();
            throw new Exception("Error to annotate tree: " + e.getMessage() + "\nPlease check the tree log file format.");
        }
        
        
        
        // Concerted evolution events
        processMetaData(targetTree.getRoot());
        
        
        // Print the tree to output
        try {
            final PrintStream stream = new PrintStream(new FileOutputStream(outputInput.get()));
            targetTree.init(stream);
            stream.println();
            stream.print("tree TREE1 = ");
            int[] dummy = new int[1];
            String newick = targetTree.getRoot().toSortedNewick(dummy, true);
            stream.print(newick);
            stream.println(";");
            targetTree.close(stream);
            stream.println();
        } catch (Exception e) {
            Log.err.println("Error to write annotated tree file: " + e.getMessage());
            return;
        }
		
        
        
		
	}
	
	/**
	 * Put metadata onto node
	 * @param node
	 */
	 private void processMetaData(Node node) {
		for (Node child : node.getChildren()) {
			processMetaData(child);
		}
		Set<String> metaDataNames = node.getMetaDataNames(); 
		if (metaDataNames != null && !metaDataNames.isEmpty()) {
			String metadata = "";
			for (String name : metaDataNames) {
				Object value = node.getMetaData(name);
				metadata += name + "=";
				if (value instanceof Object[]) {
					Object [] values = (Object[]) value;
					metadata += "{";
					for (int i = 0; i < values.length; i++) {
						metadata += values[i].toString();
						if (i < values.length - 1) {
							metadata += ",";
						}
					}
					metadata += "}";
				} else {
					 metadata += value.toString();
				}
				metadata += ",";
			}
			metadata = metadata.substring(0, metadata.length() - 1);
			node.metaDataString = metadata;
		}		
	}


	

	/**
	 * Get MCC tree (borrowed from TreeAnnotator)
	 * @param cladeSystem
	 * @param useSumCladeCredibility
	 * @return
	 * @throws IOException
	 */
    protected Tree summarizeTrees(PrunedTreeSet treeSet, CladeSystem cladeSystem, List<Boolean> include) throws IOException  {

        Tree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        System.out.println("Analyzing " + totalTreesUsed + " trees...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");

        int stepSize = Math.max(totalTreesUsed / 60, 1);
        int reported = 0;

        int counter = 0;
        treeSet.reset();
        while (treeSet.hasNext()) {
        		Tree tree = treeSet.next();
        		if (!include.get(counter)) {
	        		counter++;
	        		continue;
        		}
            	double score = cladeSystem.getLogCladeCredibility(tree.getRoot(), null); //
			  if (score > bestScore) {
			      bestTree = tree;
			      bestScore = score;
			  }
			  if (counter % stepSize == 0 && reported < 61) {
				  while (1000*reported < 61000 * (counter + 1) / totalTreesUsed) {
					  System.out.print("*");
			          reported++;
				  }
				  System.out.flush();
			  }
			  counter++;
        }
        System.out.println();
        System.out.println();
    	System.out.println("Highest Log Clade Credibility: " + bestScore);

        return bestTree;
    }
	
    
	
	
    /**
     * Find the most commonly occurring clade set and remove all other occurrences
     * @param treeSet
     * @return
     * @throws IOException 
     */
    private List<Boolean> filterTreeset(PrunedTreeSet treeSet) throws IOException {
    	
    	
    	Map<String, Integer> cladeCounts = new HashMap<>();
    	List<String> clades = new ArrayList<>();
    	
    			
    	treeSet.reset();
    	while (treeSet.hasNext()) {
    		Tree tree = treeSet.next();
            tree.getLeafNodeCount();
            
            
            // Get clade string for this tree
            List<Node> tips = tree.getRoot().getAllLeafNodes();
    		if (tips.isEmpty()) tips.add(tree.getRoot());
            List<String> tipsClade = new ArrayList<>();
    		for (int tipNr = 0; tipNr < tips.size(); tipNr++) {
    			tipsClade.add(tips.get(tipNr).getID());
    		}
    		Collections.sort(tipsClade);
    		String clade = String.join(" + ", tipsClade);
    		
    		
    		// Increment count
    		int count = 1;
    		if (cladeCounts.containsKey(clade)) {
    			count += cladeCounts.get(clade);
    		}
    		cladeCounts.put(clade, count);
    		clades.add(clade);
    		
    		
        } 
    	
    	
    	// Which clade is the most popular?
    	int maxCount = 0;
    	String mostPopularClade = null;
    	for (String clade : cladeCounts.keySet()) {
    		if (cladeCounts.get(clade) > maxCount) {
    			mostPopularClade = clade;
    			maxCount = cladeCounts.get(clade);
    		}
    	}
    	
    	Log.warning("Identified most commonly occurring clade in " + maxCount + " trees");
    	
    	
    	// Filter
    	int n = 0;
    	List<Boolean> include = new ArrayList<>();
    	for (String clade : clades) {
    		boolean includeTree = clade.equals(mostPopularClade);
    		include.add(includeTree);
    		if (includeTree) n++;
    	}
    	
    	//Log.warning(n + " trees");
    	
    	
		// TODO Auto-generated method stub
		return include;
	}

	/**
     * Annotate a tree using attributes (borrowed from TreeAnnotator)
     * @param tree
     */
	private void setupAttributes(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            Set<String> iter = node.getMetaDataNames();
            if (iter != null) {
            	for (String name : iter) {
                    attributeNames.add(name);
                }
            }
        }

    }

	
	
	

    protected void annotateTree(CladeSystem cladeSystem, Node node, BitSet bits) {

        BitSet bits2 = new BitSet();

        if (node.isLeaf()) {

            int index = cladeSystem.getTaxonIndex(node);
            bits2.set(2*index);

            annotateNode(cladeSystem, node, bits2, true);
        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                annotateTree(cladeSystem, node1, bits2);
            }

            for (int i=1; i<bits2.length(); i=i+2) {
                bits2.set(i, false);
            }
           

            annotateNode(cladeSystem, node, bits2, false);
        }

        if (bits != null) {
            bits.or(bits2);
        }
    }

    
    
    /**
     * Annotate a tree (borrowed from TreeAnnotator)
     * @param cladeSystem
     * @param node
     * @param bits
     * @param isTip
     * @param heightsOption
     */
    protected void annotateNode(CladeSystem cladeSystem, Node node, BitSet bits, boolean isTip) {
        CladeSystem.Clade clade = cladeSystem.getCladeMap().get(bits);
        assert clade != null : "Clade missing?";

        boolean filter = false;
        if (!isTip) {
            final double posterior = clade.getCredibility();
            node.setMetaData("posterior", posterior);
        }
        
        
        // Find all concerted shift events for this clade
        List<String> taxa = new ArrayList<>();
        if (node.isLeaf()) taxa.add(node.getID());
        else for (Node leaf : node.getAllLeafNodes()) {
        	taxa.add(leaf.getID());
        }
        Collections.sort(taxa);
        String cladeStr = String.join(" + ", taxa);
        int eventNr = 1;
        //Log.warning("Clade " + cladeStr);
        
        


        int i = 0;
        for (String attributeName : attributeNames) {
        	
        	

            if (clade.getAttributeValues() != null && clade.getAttributeValues().size() > 0) {
                double[] values = new double[clade.getAttributeValues().size()];

                HashMap<Object, Integer> hashMap = new HashMap<>();

                Object[] v = clade.getAttributeValues().get(0);
                if (v[i] != null) {

                    final boolean isHeight = attributeName.equals("height");
                    boolean isBoolean = v[i] instanceof Boolean;

                    boolean isDiscrete = v[i] instanceof String;

                    if (forceIntegerToDiscrete && v[i] instanceof Integer) isDiscrete = true;

                    double minValue = Double.MAX_VALUE;
                    double maxValue = -Double.MAX_VALUE;

                    final boolean isArray = v[i] instanceof Object[];
                    boolean isDoubleArray = isArray && ((Object[]) v[i])[0] instanceof Double;
                    // This is Java, friends - first value type does not imply all.
                    if (isDoubleArray) {
                        for (Object n : (Object[]) v[i]) {
                            if (!(n instanceof Double)) {
                                isDoubleArray = false;
                                break;
                            }
                        }
                    }
                    // todo Handle other types of arrays

                    double[][] valuesArray = null;
                    double[] minValueArray = null;
                    double[] maxValueArray = null;
                    int lenArray = 0;

                    if (isDoubleArray) {
                        lenArray = ((Object[]) v[i]).length;

                        valuesArray = new double[lenArray][clade.getAttributeValues().size()];
                        minValueArray = new double[lenArray];
                        maxValueArray = new double[lenArray];

                        for (int k = 0; k < lenArray; k++) {
                            minValueArray[k] = Double.MAX_VALUE;
                            maxValueArray[k] = -Double.MAX_VALUE;
                        }
                    }

                    for (int j = 0; j < clade.getAttributeValues().size(); j++) {
                        Object value = clade.getAttributeValues().get(j)[i];
                        if (isDiscrete) {
                            final Object s = value;
                            if (s == null) continue;
                            if (hashMap.containsKey(s)) {
                                hashMap.put(s, hashMap.get(s) + 1);
                            } else {
                                hashMap.put(s, 1);
                            }
                        } else if (isBoolean) {
                            values[j] = (((Boolean) value) ? 1.0 : 0.0);
                        } else if (isDoubleArray) {
                            // Forcing to Double[] causes a cast exception. MAS
                            try {
                                Object[] array = (Object[]) value;
                                for (int k = 0; k < lenArray; k++) {
                                    valuesArray[k][j] = ((Double) array[k]);
                                    if (valuesArray[k][j] < minValueArray[k]) minValueArray[k] = valuesArray[k][j];
                                    if (valuesArray[k][j] > maxValueArray[k]) maxValueArray[k] = valuesArray[k][j];
                                }
                            } catch (Exception e) {
                                // ignore
                            }
                        } else {
                            // Ignore other (unknown) types
                            if (value instanceof Number) {
                                values[j] = ((Number) value).doubleValue();
                                if (values[j] < minValue) minValue = values[j];
                                if (values[j] > maxValue) maxValue = values[j];
                            }
                        }
                    }
                   

                    if (!filter) {
                        boolean processed = false;
                       

                        if (!processed) {
                            if (!isDiscrete) {
                                if (!isDoubleArray) {
                                	
                                	String attrNameLabel = attributeName;
                                    annotateMeanAttribute(node, attrNameLabel, values);
                                } else {
                                    for (int k = 0; k < lenArray; k++) {
                                        annotateMeanAttribute(node, attributeName + (k + 1), valuesArray[k]);
                                    }
                                }
                            } else {
                                annotateModeAttribute(node, attributeName, hashMap);
                                annotateFrequencyAttribute(node, attributeName, hashMap);
                            }
                            if (!isBoolean && minValue < maxValue && !isDiscrete && !isDoubleArray) {
                                // Basically, if it is a boolean (0, 1) then we don't need the distribution information
                                // Likewise if it doesn't vary.
                                annotateMedianAttribute(node, attributeName + "_median", values);
                                annotateHPDAttribute(node, attributeName + "_95%_HPD", 0.95, values);
                                annotateRangeAttribute(node, attributeName + "_range", values);
                            }

                            if (isDoubleArray) {
                                String name = attributeName;
                                // todo
//                                    if (name.equals(location1Attribute)) {
//                                        name = locationOutputAttribute;
//                                    }
                                boolean want2d = lenArray == 2;
                                if (name.equals("dmv")) {  // terrible hack
                                    want2d = false;
                                }
                                for (int k = 0; k < lenArray; k++) {
                                    if (minValueArray[k] < maxValueArray[k]) {
                                        annotateMedianAttribute(node, name + (k + 1) + "_median", valuesArray[k]);
                                        annotateRangeAttribute(node, name + (k + 1) + "_range", valuesArray[k]);
                                        if (!want2d)
                                            annotateHPDAttribute(node, name + (k + 1) + "_95%_HPD", 0.95, valuesArray[k]);
                                    }
                                }
                                // 2D contours
                                if (want2d) {

                                    boolean variationInFirst = (minValueArray[0] < maxValueArray[0]);
                                    boolean variationInSecond = (minValueArray[1] < maxValueArray[1]);

                                    if (variationInFirst && !variationInSecond)
                                        annotateHPDAttribute(node, name + "1" + "_95%_HPD", 0.95, valuesArray[0]);

                                    if (variationInSecond && !variationInFirst)
                                        annotateHPDAttribute(node, name + "2" + "_95%_HPD", 0.95, valuesArray[1]);

                                    if (variationInFirst && variationInSecond)
                                        annotate2DHPDAttribute(node, name, "_" + (int) (100 * hpd2D) + "%HPD", hpd2D, valuesArray);
                                }
                            }
                        }
                    }
                }
            }
            i++;
        }
    }
	
	

    protected void annotateMeanAttribute(Node node, String label, double[] values) {
        double mean = DiscreteStatistics.mean(values);
        node.setMetaData(label, mean);
    }

    protected void annotateMedianAttribute(Node node, String label, double[] values) {
        double median = DiscreteStatistics.median(values);
        node.setMetaData(label, median);

    }

    protected void annotateModeAttribute(Node node, String label, HashMap<Object, Integer> values) {
        Object mode = null;
        int maxCount = 0;
        int totalCount = 0;
        int countInMode = 1;

        for (Object key : values.keySet()) {
            int thisCount = values.get(key);
            if (thisCount == maxCount) {
                // I hope this is the intention
                mode = mode.toString().concat("+" + key);
                countInMode++;
            } else if (thisCount > maxCount) {
                mode = key;
                maxCount = thisCount;
                countInMode = 1;
            }
            totalCount += thisCount;
        }
        double freq = (double) maxCount / (double) totalCount * countInMode;
        node.setMetaData(label, mode);
        node.setMetaData(label + ".prob", freq);
    }

    protected void annotateFrequencyAttribute(Node node, String label, HashMap<Object, Integer> values) {
        double totalCount = 0;
        Set<?> keySet = values.keySet();
        int length = keySet.size();
        String[] name = new String[length];
        Double[] freq = new Double[length];
        int index = 0;
        for (Object key : values.keySet()) {
            name[index] = key.toString();
            freq[index] = new Double(values.get(key));
            totalCount += freq[index];
            index++;
        }
        for (int i = 0; i < length; i++)
            freq[i] /= totalCount;

        node.setMetaData(label + ".set", name);
        node.setMetaData(label + ".set.prob", freq);
    }

    protected void annotateRangeAttribute(Node node, String label, double[] values) {
        double min = DiscreteStatistics.min(values);
        double max = DiscreteStatistics.max(values);
        node.setMetaData(label, new Object[]{min, max});
    }

    protected void annotateHPDAttribute(Node node, String label, double hpd, double[] values) {
        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(hpd * values.length);
        for (int i = 0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex + diff - 1]];
        node.setMetaData(label, new Object[]{lower, upper});
    }
	
    

    public static final String CORDINATE = "cordinates";
    protected String formattedLocation(double x) {
        return String.format("%5.2f", x);
    }

    
    protected void annotate2DHPDAttribute(Node node, String preLabel, String postLabel, double hpd, double[][] values) {
        

    	//  KernelDensityEstimator2D kde = new KernelDensityEstimator2D(values[0], values[1], N);
        //ContourMaker kde = new ContourWithSynder(values[0], values[1], N);
        boolean bandwidthLimit = false;

        ContourMaker kde = new ContourWithSynder(values[0], values[1], bandwidthLimit);

        ContourPath[] paths = kde.getContourPaths(hpd);

        node.setMetaData(preLabel + postLabel + "_modality", paths.length);

        if (paths.length > 1) {
            Log.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
            Log.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
        }

        StringBuffer output = new StringBuffer();
        int i = 0;
        for (ContourPath p : paths) {
            output.append("\n<" + CORDINATE + ">\n");
            double[] xList = p.getAllX();
            double[] yList = p.getAllY();
            StringBuffer xString = new StringBuffer("{");
            StringBuffer yString = new StringBuffer("{");
            for (int k = 0; k < xList.length; k++) {
                xString.append(formattedLocation(xList[k])).append(",");
                yString.append(formattedLocation(yList[k])).append(",");
            }
            xString.append(formattedLocation(xList[0])).append("}");
            yString.append(formattedLocation(yList[0])).append("}");

            node.setMetaData(preLabel + "1" + postLabel + "_" + (i + 1), xString);
            node.setMetaData(preLabel + "2" + postLabel + "_" + (i + 1), yString);
            i++;

        }
        
            
            
    }
    
    /**
     * Borrowed from CladeSystem
     * @param node
     * @param codes
     * @return
     */
    protected int getTreeCladeCodes(Node node, BitSet[] codes) {
        final int inode = node.getNr();
        codes[inode].clear();
        if (node.isLeaf()) {
            int index = node.getNr();
            codes[inode].set(index);
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                final Node child = node.getChild(i);
                final int childIndex = getTreeCladeCodes(child, codes);

                codes[inode].or(codes[childIndex]);
            }
        }
        return inode;
    }
	
	

    protected boolean setTreeHeightsByCA(PrunedTreeSet treeSet, Tree targetTree, List<Boolean> include) throws IOException  {
    	
        System.out.println("Setting node heights...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");

        int reportStepSize = totalTreesUsed / 60;
        if (reportStepSize < 1) reportStepSize = 1;
        int reported = 0;


        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        CladeSystem cladeSystem = new CladeSystem(targetTree);
        final int clades = cladeSystem.getCladeMap().size();

        // allocate posterior tree nodes order once
        int[] postOrderList = new int[clades];
        BitSet[] ctarget = new BitSet[clades];
        BitSet[] ctree = new BitSet[clades];

        for (int k = 0; k < clades; ++k) {
            ctarget[k] = new BitSet();
            ctree[k] = new BitSet();
        }

        getTreeCladeCodes(targetTree.getRoot(), ctarget);

        // temp collecting heights inside loop allocated once
        double[][] hs = new double[clades][totalTreesUsed]; //[treeSet.totalTrees - treeSet.burninCount];

        // heights total sum from posterior trees
        double[] ths = new double[clades];

        int totalTreesUsed = 0;

        int counter = 0;
        int treeNum = 0;
        treeSet.reset();
        while (treeSet.hasNext()) {
        	Tree tree = treeSet.next();
        	
        	if (include.get(counter)) {
        	
	            TreeUtils.preOrderTraversalList(tree, postOrderList);
	            getTreeCladeCodes(tree.getRoot(), ctree);
	            for (int k = 0; k < clades; ++k) {
	                int j = postOrderList[k];
	                for (int i = 0; i < clades; ++i) {
	                    if( CollectionUtils.isSubSet(ctarget[i], ctree[j]) ) {
	                        hs[i][treeNum] = tree.getNode(j).getHeight();
	                    }
	                }
	            }
	            for (int k = 0; k < clades; ++k) {
	                ths[k] += hs[k][treeNum];
	            }
	            totalTreesUsed += 1;
	            if (treeNum > 0 && treeNum % reportStepSize == 0 && reported < 61) {
					while (1000 * reported < 61000 * (treeNum + 1)/ this.totalTreesUsed) {
					    System.out.print("*");
					    reported++;
					}
					 System.out.flush();
	            }
            
	            treeNum++;
        	}
        	
            counter++;

        }

        targetTree.initAndValidate();

        cladeSystem.removeClades(targetTree.getRoot(), true);
        for (int k = 0; k < clades; ++k) {
            ths[k] /= totalTreesUsed;
            final Node node = targetTree.getNode(k);
            node.setHeight(ths[k]);
            String attributeName = "CAheight";
            double [] values = hs[k];
            double min = values[0];
            double max = values[0];
            for (double d : values) {
            	min = Math.min(d, min);
            	max = Math.max(d, max);
            }
            if (Math.abs(min - max) > 1e-10) {
	            annotateMeanAttribute(node, attributeName + "_mean", values);
	            annotateMedianAttribute(node, attributeName + "_median", values);
	            annotateHPDAttribute(node, attributeName + "_95%_HPD", 0.95, values);
	            annotateRangeAttribute(node, attributeName + "_range", values);
            }
        }

        assert (totalTreesUsed == this.totalTreesUsed);
        this.totalTreesUsed = totalTreesUsed;
        System.out.println();
        System.out.println();

        return true;
    }

    
    
	
	public static void main(String [] args) throws Exception {
		new Application(new PrunedTreeAnnotator(), "Pruned tree annotator", args);
	} 


}
