package uniandes.algobc.metabolites;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class Graph<K,E> {
	Map<K, List<Node<K,E>>> G = null;
	public Graph() {
		G = new HashMap<>();
	}
	public boolean addEdge(K v1, K v2, E edgeLabel) {
		if(!G.containsKey(v1)) {
			G.put(v1, new LinkedList<>());
		}
		G.get(v1).add(new Node<>(v1, v2, edgeLabel));
		return true;
	}
	public String toString() {
		String data = "";
		for(K key : G.keySet()) {
			data += G.get(key).toString().replaceAll(",", "").replaceAll("\\[","").replace("]", "");
		}
		return data;
	}
}
