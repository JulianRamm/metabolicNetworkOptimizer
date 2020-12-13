package uniandes.algobc.metabolites;

/**
 * This represents the data related to a metabolite within a reaction
 * @author Jorge Duitama
 *
 */
public class ReactionComponent {
	private Metabolite metabolite;
	private double stoichiometry;
	/**
	 * Creates a reaction component with the given data
	 * @param metabolite that participates in the reaction
	 * @param stoichiometry coefficient of the metabolite in the reaction
	 */
	public ReactionComponent(Metabolite metabolite, double stoichiometry) {
		super();
		this.metabolite = metabolite;
		this.stoichiometry = stoichiometry;
	}
	/**
	 * @return Metabolite that participates in the reaction
	 */
	public Metabolite getMetabolite() {
		return metabolite;
	}
	/**
	 * 
	 * @return steichiometry coefficient of the metabolite within the reaction
	 */
	public double getStoichiometry() {
		return stoichiometry;
	}
	
	
}
