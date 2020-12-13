package uniandes.algobc.metabolites;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
/**
 * Represents a metabolic network of reactions on metabolites
 * @author Jorge Duitama
 */
public class MetabolicNetwork {

	private Map<String,Enzyme> enzymes = new TreeMap<String,Enzyme>(); 
	private Map<String,Metabolite> metabolites = new TreeMap<String,Metabolite>();
	private Set<String> compartments = new TreeSet<String>();
	private Map<String,Reaction> reactions = new TreeMap<String,Reaction>();
	private Graph<Metabolite, Integer> metaboliteGraph = new Graph<Metabolite, Integer>();
	private int[][] adjacencyMatrix;
	
	/**
	 * Adds a new gene product that can catalyze reactions
	 * @param product New gene product
	 */
	public void addEnzyme(Enzyme enzyme) {
		enzymes.put(enzyme.getId(), enzyme);
	}
	/**
	 * Adds a new metabolite. If a metabolite with the given name is already added, it 
	 * @param metabolite New metabolite
	 */
	public void addMetabolite(Metabolite metabolite) {
		metabolites.put(metabolite.getId(), metabolite);
		compartments.add(metabolite.getCompartment());
	}
	/**
	 * Adds a new reaction
	 * @param r New reaction between metabolites
	 */
	public void addReaction(Reaction r) {
		reactions.put(r.getId(),r);
	}
	
	public void deleteReaction(Reaction r) {
		reactions.remove(r.getId());
	}
	public void deleteMetabolite(Metabolite m) {
		metabolites.remove(m.getId());
	}
	/**
	 * Returns the gene product with the given id
	 * @param id of the product to search
	 * @return GeneProduct with the given id
	 */
	public Enzyme getEnzyme (String id) {
		return enzymes.get(id);
	}
	/**
	 * Returns the metabolite with the given id
	 * @param id of the metabolite to search
	 * @return Metabolite with the given id
	 */
	public Metabolite getMetabolite (String id) {
		return metabolites.get(id);
	}
	/**
	 * @return List of metabolites in the network
	 */
	public List<Metabolite> getMetabolitesList() {
		return new ArrayList<Metabolite>(metabolites.values());
	}
	/**
	 * @return List of reactions in the network
	 */
	public List<Reaction> getReactionsList () {
		return new ArrayList<Reaction>(reactions.values());
	}
	public ArrayList<Metabolite> getReactantsMetabolites(){
		ArrayList<Metabolite> result = new ArrayList<Metabolite>();
		for (Map.Entry<String,Reaction> entry : reactions.entrySet()) { 
				for(ReactionComponent reactant: entry.getValue().getReactants()) {
					result.add(reactant.getMetabolite());
				}
    	} 
		return result;
	}
	public ArrayList<Metabolite> getProductsMetabolites(){
		ArrayList<Metabolite> result = new ArrayList<Metabolite>();
		for (Map.Entry<String,Reaction> entry : reactions.entrySet()) { 
				for(ReactionComponent product: entry.getValue().getProducts()) {
					result.add(product.getMetabolite());
				}
    	} 
		return result;
	}
	public void createMetaboliteGraph() {
		ArrayList<Metabolite> reactants = this.getReactantsMetabolites();
		ArrayList<Metabolite> products = this.getProductsMetabolites();
		for( Metabolite m: reactants) {
			for(Metabolite m1: products) {
				int weight = findNumberReactions(m, m1);
				if(weight > 1 ) {
					metaboliteGraph.addEdge(m, m1, weight);
				}
			}
		}
	}
	public int[][] createAdjacencyMatrix() {
		ArrayList<Metabolite> reactants = this.getReactantsMetabolites();
		ArrayList<Metabolite> products = this.getProductsMetabolites();
		adjacencyMatrix = new int[reactants.size()][products.size()];
		for( int i = 0; i<reactants.size(); i++ ) {
			for(int j = 0; j<products.size(); j++) {
				adjacencyMatrix[i][j] = findNumberReactions(reactants.get(i), products.get(j));
			}
		}
		return adjacencyMatrix;
	}

	public int findNumberReactions(Metabolite reactant, Metabolite product) {
		int res = 0;
		for (Map.Entry<String,Reaction> entry : reactions.entrySet()) {
			for(ReactionComponent rc: entry.getValue().getReactants()) {
				if(rc.getMetabolite().equals(reactant)) {
					for(ReactionComponent rc1: entry.getValue().getProducts()) {
						if(rc1.getMetabolite().equals(product)) {
							res++;
						}
					}
				}
			}
		}

		return res;
	}
	public MetabolicNetwork optimize(List<String> metaboliteIds, List<String> reactionIds, List<Inequality> functions, int mDof, int nMin) {
		Set<Metabolite> protectedMetabolites = new HashSet<Metabolite>();
		Set<Reaction> protectedReactions = new HashSet<Reaction>();
		Map<Inequality, Boolean> protectedFunctionsAndPhenotypes = new HashMap<Inequality, Boolean>();
		for (Reaction r: getReactionsList()) {
			if(reactionIds.contains(r.getId())) {
				protectedReactions.add(r);
			}
		}
		for (Metabolite m: getMetabolitesList()) {
			if(metaboliteIds.contains(m.getId())) {
				protectedMetabolites.add(m);
			}
		}
		for (Inequality ineq: functions) {
			protectedFunctionsAndPhenotypes.put(ineq, false);
		}
		return new MetabolicNetworkOptimizer(this, protectedMetabolites, protectedReactions, protectedFunctionsAndPhenotypes, mDof, nMin).optimize();
	}
	public String getData() {
		return metaboliteGraph.toString();
	}
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		System.out.println("Enzymes");
		for(Enzyme enzyme:network.enzymes.values()) {
			System.out.println(enzyme.getId()+" "+enzyme.getName());
		}
		System.out.println();
		
		List<Metabolite> metabolitesList = network.getMetabolitesList();
		System.out.println("Loaded "+metabolitesList.size()+" metabolites: ");
		for(Metabolite m:metabolitesList) {
			System.out.println(m.getId()+" "+m.getName()+" "+m.getCompartment()+" "+m.getChemicalFormula());
		}
		System.out.println();
		List<Reaction> reactions = network.getReactionsList();
		System.out.println("Loaded "+reactions.size()+" reactions");
		for(Reaction r:reactions) {
			System.out.println(r.getId()+" "+r.getName()+" "+r.getReactants().size()+" "+r.getProducts().size()+" "+r.getEnzymes().size()+" "+r.getLowerBoundFlux()+" "+r.getUpperBoundFlux());
		}
		//Ejemplo de lista de metabolitos: { "Pkm2", "sada", "asdas"};
		//Ejemplo de lista de reacciones: {"asdjs", "reactj", "asd"};
		int mDof = 100; // mínimos grados de libertad a utilizar
		int nMin = 100; // mínimo númeor de reacciones
		List<String> protectedMetabolites = Arrays.asList(new String[] { "Pkm2", "sada", "asdas"}); //Lista de los nombres de los metabolitos protegidos
		List<String> protectedReactions = Arrays.asList(new String[] {"asdjs", "reactj", "asd"}); // Lista de los nombres las reacciones protegidas
		List<Inequality> protectedFunctionsAndPhenotypes = new ArrayList<Inequality>();
		Inequality ejemplo = new Inequality();
		ejemplo.createInequality(-1.0, -0.999, 0.5, -1000.0, 1000.0, -1000.0, 1000.0); //Ejmplo de cómo crear una desigualdad
		protectedFunctionsAndPhenotypes.add(ejemplo);
		MetabolicNetwork optimizedNetwork = network.optimize(protectedMetabolites, protectedReactions, protectedFunctionsAndPhenotypes, mDof, nMin); //Red metabólica reducida
		optimizedNetwork.createMetaboliteGraph();
		try {
		      FileWriter myWriter = new FileWriter("data/grafo.txt");
		      myWriter.write(optimizedNetwork.getData());
		      
		      myWriter.close();
		      System.out.println("Successfully wrote to the file.");
		    } catch (IOException e) {
		      System.out.println("An error occurred.");
		      e.printStackTrace();
		    }
		System.out.println(optimizedNetwork.metaboliteGraph);
		System.out.println(optimizedNetwork.getReactionsList().size());
		System.out.println(optimizedNetwork.getMetabolitesList().size());
	}
}
