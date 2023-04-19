package indeldollo.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;

/**
 * Code adapted from MemoryFriendlyTreeSet
 * @author jdou557
 *
 */


@Description("Like TreeSet, but more robust at tree parsing and does not require all trees to have the same taxonset")
public class PrunedTreeSet {
	
	static PrintStream progressStream = Log.err;
    String inputFileName;
    int burninCount = 0;
    int totalTrees = 0;
    boolean isNexus = true;
    
    
    int current = 0;
	int lineNr;
    public Map<String, String> translationMap = null;
    public List<String> taxa;

    // label count origin for NEXUS trees
    int origin = -1;

    BufferedReader fin;
    
    
    public PrunedTreeSet(String inputFileName, int burninPercentage) throws IOException  {
		this.inputFileName = inputFileName;
		countTrees(burninPercentage);
        fin = new BufferedReader(new FileReader(inputFileName));
	}


    /** determine number of trees in the file,
	 * and number of trees to skip as burnin
	 * @throws IOException
	 * @throws FileNotFoundException **/
	void countTrees(int burninPercentage) throws IOException  {
        BufferedReader fin = new BufferedReader(new FileReader(new File(inputFileName)));
        if (!fin.ready()) {
        	throw new IOException("File appears empty");
        }
    	String str = fin.readLine();
        if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
        	// the file contains a list of Newick trees instead of a list in Nexus format
        	isNexus = false;
        	if (str.trim().length() > 0) {
        		totalTrees = 1;
        	}
        }
        while (fin.ready()) {
        	str = fin.readLine();
            if (isNexus) {
                if (str.trim().toLowerCase().startsWith("tree ")) {
                	totalTrees++;
                }
            } else if (str.trim().length() > 0) {
        		totalTrees++;
            }
        }
        fin.close();

        burninCount = Math.max(0, (burninPercentage * totalTrees)/100);

        progressStream.println("Processing " + (totalTrees - burninCount) + " trees from file" +
                (burninPercentage > 0 ? " after ignoring first " + burninPercentage + "% = " + burninCount + " trees." : "."));
	}

	public boolean hasNext() {
		return current < totalTrees;
	}

	public Tree next() throws IOException {
		String str = nextLine();
		if (!isNexus) {
            PrunedTreeParser treeParser;

            if (origin != -1) {
                treeParser = new PrunedTreeParser(taxa, str, origin, false);
            } else {
                try {
                    treeParser = new PrunedTreeParser(taxa, str, 0, false);
                } catch (ArrayIndexOutOfBoundsException e) {
                    treeParser = new PrunedTreeParser(taxa, str, 1, false);
                }
            }
            return treeParser;
		}
		
        // read trees from NEXUS file
        if (str.trim().toLowerCase().startsWith("tree ")) {
        	current++;
            final int i = str.indexOf('(');
            if (i > 0) {
                str = str.substring(i);
            }
            PrunedTreeParser treeParser;

            if (origin != -1) {
                treeParser = new PrunedTreeParser(taxa, str, origin, false);
            } else {
                try {
                    treeParser = new PrunedTreeParser(taxa, str, 0, false);
                } catch (ArrayIndexOutOfBoundsException e) {
                    treeParser = new PrunedTreeParser(taxa, str, 1, false);
                }
            }

            if (translationMap != null) treeParser.translateLeafIds(translationMap);

            return treeParser;
        }
		return null;
	}
	
	/**
     * read next line from Nexus file that is not a comment and not empty 
     * @throws IOException *
     */
    String nextLine() throws IOException  {
        String str = readLine();
        if (str == null) {
            return null;
        }
        if (str.matches("^\\s*\\[.*")) {
            final int start = str.indexOf('[');
            int end = str.indexOf(']', start);
            while (end < 0) {
                str += readLine();
                end = str.indexOf(']', start);
            }
            str = str.substring(0, start) + str.substring(end + 1);
            if (str.matches("^\\s*$")) {
                return nextLine();
            }
        }
        if (str.matches("^\\s*$")) {
            return nextLine();
        }
        return str;
    }

    /**
     * read line from nexus file *
     */
    String readLine() throws IOException {
        if (!fin.ready()) {
            return null;
        }
        lineNr++;
        return fin.readLine();
    }

	

	public void reset() throws IOException {
		current = 0;
        fin = new BufferedReader(new FileReader(new File(inputFileName)));
        lineNr = 0;
        try {
            while (fin.ready()) {
                final String str = nextLine();
                if (str == null) {
                    return;
                }
                final String lower = str.toLowerCase();
                if (lower.matches("^\\s*begin\\s+trees;\\s*$")) {
                    parseTreesBlock();
                    return;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Around line " + lineNr + "\n" + e.getMessage());
        }
		
	}
	
	

    private void parseTreesBlock() throws IOException  {
        // read to first non-empty line within trees block
    	fin.mark(1024*1024);
    	int lineNr = this.lineNr;
        String str = readLine().trim();
        while (str.equals("")) {
        	fin.mark(1024*1024);
        	lineNr = this.lineNr;
            str = readLine().trim();
        }

        // if first non-empty line is "translate" then parse translate block
        if (str.toLowerCase().contains("translate")) {
            translationMap = parseTranslateBlock();
            origin = getIndexedTranslationMapOrigin(translationMap);
            if (origin != -1) {
                taxa = getIndexedTranslationMap(translationMap, origin);
            }
        } else {
        	this.lineNr = lineNr;
        	fin.reset();
        }
        // we got to the end of the translate block
        // read burninCount trees
        current = 0;
        while (current < burninCount && fin.ready()) {
			str = nextLine();
            if (str.trim().toLowerCase().startsWith("tree ")) {
            	current++;
            }
        }
    }

    private List<String> getIndexedTranslationMap(final Map<String, String> translationMap, final int origin) {

        //System.out.println("translation map size = " + translationMap.size());

        final String[] taxa = new String[translationMap.size()];

        for (final String key : translationMap.keySet()) {
            taxa[Integer.parseInt(key) - origin] = translationMap.get(key);
        }
        return Arrays.asList(taxa);
    }

    /**
     * @param translationMap
     * @return minimum key value if keys are a contiguous set of integers starting from zero or one, -1 otherwise
     */
    private int getIndexedTranslationMapOrigin(final Map<String, String> translationMap) {

        final SortedSet<Integer> indices = new java.util.TreeSet<>();

        int count = 0;
        for (final String key : translationMap.keySet()) {
            final int index = Integer.parseInt(key);
            indices.add(index);
            count += 1;
        }
        if ((indices.last() - indices.first() == count - 1) && (indices.first() == 0 || indices.first() == 1)) {
            return indices.first();
        }
        return -1;
    }

    /**
     * @return a map of taxa translations, keys are generally integer node number starting from 1
     *         whereas values are generally descriptive strings.
     * @throws IOException
     */
    private Map<String, String> parseTranslateBlock() throws IOException {

        final Map<String, String> translationMap = new HashMap<>();

        String line = readLine();
        final StringBuilder translateBlock = new StringBuilder();
        while (line != null && !line.trim().toLowerCase().equals(";")) {
            translateBlock.append(line.trim());
            line = readLine();
        }
        final String[] taxaTranslations = translateBlock.toString().split(",");
        for (final String taxaTranslation : taxaTranslations) {
            final String[] translation = taxaTranslation.split("[\t ]+");
            if (translation.length == 2) {
                translationMap.put(translation[0], translation[1]);
//                System.out.println(translation[0] + " -> " + translation[1]);
            } else {
                Log.err.println("Ignoring translation:" + Arrays.toString(translation));
            }
        }
        return translationMap;
    }

	

}
