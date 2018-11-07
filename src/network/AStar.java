package network;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.planargraph.DirectedEdgeStar;
import com.vividsolutions.jts.planargraph.Node;

import sim.field.network.Edge;
import sim.field.network.Network;
import sim.util.geo.GeomPlanarGraphDirectedEdge;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * AStar.java
 *
 * Copyright 2011 by Sarah Wise, Mark Coletti, Andrew Crooks, and George Mason
 * University.
 *
 * Licensed under the Academic Free License version 3.0
 *
 * See the file "LICENSE" for more information
 *
 * $Id: AStar.java 842 2012-12-18 01:09:18Z mcoletti $
 */
public class AStar {

	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// PARAMETERS ///////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	public ArrayList<Edge> astarWeightedPath(GeoNode start, GeoNode goal, Network network, HashSet<Object> weighted,
			double weight) {
		// initial check
		if (start == null || goal == null) {
			System.out.println("Error: invalid node provided to AStar");
			return null;
		}

		// if they're the same place, the path is empty but certainly exists
		if (start == goal)
			return new ArrayList<Edge>();

		// set up the containers for the result
		// ArrayList<GeomPlanarGraphDirectedEdge> result =
		// new ArrayList<GeomPlanarGraphDirectedEdge>();

		// containers for the metainformation about the Nodes relative to the
		// A* search
		HashMap<GeoNode, AStarNodeWrapper> foundNodes = new HashMap<GeoNode, AStarNodeWrapper>();

		AStarNodeWrapper startNode = new AStarNodeWrapper(start);
		AStarNodeWrapper goalNode = new AStarNodeWrapper(goal);
		foundNodes.put(start, startNode);
		foundNodes.put(goal, goalNode);

		startNode.gx = 0;
		startNode.hx = heuristic(start, goal);
		startNode.fx = heuristic(start, goal);

		// A* containers: nodes to be investigated, nodes that have been investigated
		ArrayList<AStarNodeWrapper> closedSet = new ArrayList<AStarNodeWrapper>(),
				openSet = new ArrayList<AStarNodeWrapper>();
		openSet.add(startNode);

		while (openSet.size() > 0) {
			// while there are reachable nodes to investigate
			AStarNodeWrapper x = findMin(openSet);
			// find the shortest path so far
			if (x.node == goal) {
				// we have found the shortest possible path to the goal!
				// Reconstruct the path and send it back.
				return reconstructPath(goalNode);
			}
			openSet.remove(x);
			// maintain the lists
			closedSet.add(x);

			// check all the edges out from this Node
			// DirectedEdgeStar des = x.node.getOutEdges();
			for (Object o : network.getEdgesOut(x.node)) {
				Edge l = (Edge) o;
				GeoNode next = null;
				next = (GeoNode) l.getOtherNode(x.node);

				// get the A* meta information about this Node
				AStarNodeWrapper nextNode;
				if (foundNodes.containsKey(next)) {
					nextNode = foundNodes.get(next);
				} else {
					nextNode = new AStarNodeWrapper(next);
					foundNodes.put(next, nextNode);
				}

				if (closedSet.contains(nextNode)) {
					// it has already been considered
					continue;
				}

				// otherwise evaluate the cost of this node/edge combo
				double edge_factor = 1, node_factor = 1;
				if (weighted.contains(l))
					edge_factor = weight;
				if (weighted.contains(next))
					node_factor = weight;
				double tentativeCost = node_factor * (x.gx + edge_factor * length(l));
				boolean better = false;

				if (!openSet.contains(nextNode)) {
					openSet.add(nextNode);
					nextNode.hx = heuristic(next, goal);
					better = true;
				} else if (tentativeCost < nextNode.gx) {
					better = true;
				}

				// store A* information about this promising candidate node
				if (better) {
					nextNode.cameFrom = x;
					nextNode.edgeFrom = l;
					nextNode.gx = tentativeCost;
					nextNode.fx = nextNode.gx + nextNode.hx;
				}
			}
		}
		// return result;
		System.out.println("A* Problem: graph has only " + closedSet.size() + " nodes associated with it");
		return null;
	}

	/**
	 * ///////////////////////////// Path Array //////////////////////////////////
	 * Takes the information about the given node n and returns the path that found
	 * it.
	 * 
	 * @param n
	 *            the end point of the path
	 * @return an ArrayList of GeomPlanarGraphDirectedEdges that lead from the given
	 *         Node to the Node from which the serach began
	 */
	ArrayList<Edge> reconstructPath(AStarNodeWrapper n) {
		ArrayList<Edge> result = new ArrayList<Edge>();
		AStarNodeWrapper x = n;
		while (x.cameFrom != null) {
			result.add(0, x.edgeFrom); // add this edge to the front of the list
			x = x.cameFrom;
		}

		// Collections.reverse(result);
		if (result.size() < 1)
			System.out.println("Stupid path...");
		return result;
	}

	/**
	 * /////////////////////////// Euclidean Distance ////////////////////////////
	 * Measure of the estimated distance between two Nodes. Extremely basic, just
	 * Euclidean distance as implemented here.
	 * Takes into account whether either of the GeoNodes entails a delay
	 * 
	 * @param x
	 * @param y
	 * @return notional "distance" between the given nodes.
	 */
	double heuristic(GeoNode x, GeoNode y) {
		Coordinate xnode = x.geometry.getCoordinate();
		Coordinate ynode = y.geometry.getCoordinate();
		int nodeCost = 0;
		if (x.hasAttribute("delay"))
			nodeCost += x.getIntegerAttribute("delay");
		if (y.hasAttribute("delay"))
			nodeCost += y.getIntegerAttribute("delay");

		return nodeCost + Math.sqrt(Math.pow(xnode.x - ynode.x, 2) + Math.pow(xnode.y - ynode.y, 2));
	}
	
	double length(Edge e) {
		Coordinate xnode = ((GeoNode) e.from()).geometry.getCoordinate();
		Coordinate ynode = ((GeoNode) e.to()).geometry.getCoordinate();
		return Math.sqrt(Math.pow(xnode.x - ynode.x, 2) + Math.pow(xnode.y - ynode.y, 2));
	}

	/**
	 * //////////////////////// Nodes to Consider ///////////////////////////////
	 * Considers the list of Nodes open for consideration and returns the node with
	 * minimum fx value
	 * 
	 * @param set
	 *            list of open Nodes
	 * @return
	 */
	AStarNodeWrapper findMin(ArrayList<AStarNodeWrapper> set) {
		// System.out.println("Finding minimum route...");
		double min = Double.MAX_VALUE;
		AStarNodeWrapper minNode = null;
		for (AStarNodeWrapper n : set) {
			if (n.fx < min) {
				min = n.fx;
				minNode = n;
			}
		}
		return minNode;
	}

	/**
	 * 
	 * /////////////////////////// A* Node Meta Info ///////////////////////////// A
	 * wrapper to contain the A* meta information about the Nodes
	 *
	 */
	class AStarNodeWrapper {
		// the underlying Node associated with the metainformation
		GeoNode node;
		// the Node from which this Node was most profitably linked
		AStarNodeWrapper cameFrom;
		// the edge by which this Node was discovered
		Edge edgeFrom;
		double gx, hx, fx;

		public AStarNodeWrapper(GeoNode next) {
			node = next;
			gx = 0;
			hx = 0;
			fx = 0;
			cameFrom = null;
			edgeFrom = null;
		}

		public int hashCode() {
			return node.hashCode();
		}
	}
}
