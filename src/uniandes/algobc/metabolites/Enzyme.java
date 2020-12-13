package uniandes.algobc.metabolites;

/**
 * Represents an enzyme that participates in one reaction
 * @author Jorge Duitama
 */
public class Enzyme {
	private String id;
	private String name;
	private String label;
	
	/**
	 * Builds a new enzyme
	 * @param id of the enzyme
	 * @param name of the enzyme
	 */
	public Enzyme(String id, String name) {
		super();
		this.id = id;
		this.name = name;
	}
	/**
	 * @return String enzyme id
	 */
	public String getId() {
		return id;
	}
	/**
	 * @return String enzyme name
	 */
	public String getName() {
		return name;
	}
	/**
	 * @return String enzyme label
	 */
	public String getLabel() {
		return label;
	}
	/**
	 * Changes the label
	 * @param label new enzyme label
	 */
	public void setLabel(String label) {
		this.label = label;
	}
}
