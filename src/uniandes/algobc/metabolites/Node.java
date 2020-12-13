package uniandes.algobc.metabolites;

public class Node<K,E> {
	
	K key = null;
	K endNode = null;
	E edgeLabel = null;
	
	public Node(K key, K endNode, E edgeLabel) {
		this.key = key;
		this.endNode = endNode;
		this.edgeLabel = edgeLabel;
	}
	
	@Override
	public boolean equals(Object node1) {
		return key.equals(node1);
	}
	@Override
	public int hashCode() {
		return key.hashCode();
	}
	public String toString() {
		return  key + "\t"+ endNode +"\t"+edgeLabel + "\n";
	}
	
}
