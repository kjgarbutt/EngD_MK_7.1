package agents;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Random;

import network.AStar;
import network.GeoNode;
import sim.MK_7_1;
import sim.Status;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.engine.Stoppable;
import sim.field.network.Edge;
import sim.field.network.Network;
import sim.portrayal.DrawInfo2D;
import sim.util.geo.GeomPlanarGraphDirectedEdge;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import sim.util.geo.PointMoveTo;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.linearref.LengthIndexedLine;
import com.vividsolutions.jts.planargraph.Node;

import network.ListEdge;

/**
 *
 * "MK_7_1" is iteration 7.1 of my EngD project model. It varies little from
 * MK_7, and just has updated details to include Gloucestershire date. It was
 * originally adapted from the MASON demo, "Gridlock", made by Sarah Wise, Mark
 * Coletti, and Andrew Crooks.
 * 
 * The model reads a number of GIS shapefiles and displays a road network, two
 * Environment Agency flood maps and a bespoke Open Source Vulnerability Index
 * (OSVI). The model reads in a .CSV and generates a predetermined number of
 * agents with set characteristics. The agents are placed on the road network
 * and are located at a Red Cross office. The model reads a separate .CSV and
 * assigns goal locations to each agent at random from a predetermined list. The
 * agents are assigned speeds at random. Once the model is started, the agents
 * move from A to B, then they change direction and head back to their start
 * position. The process repeats until the user quits.
 *
 * @author KJGarbutt
 *
 */
public final class Agent extends TrafficAgent implements Steppable {
	private static final long serialVersionUID = -1113018274619047013L;

	////////// Objects ///////////////////////////////////////
	MK_7_1 world;

	Stoppable stopper = null;
	////////// Activities ////////////////////////////////////
	int currentActivity = 0;

	public static int activity_travel = 1;
	public static int activity_dist = 2;
	public static int activity_replen = 3;

	public boolean inbound = false;
	public boolean outbound = false;
	public boolean distributing = false;
	public boolean replenishing = false;

	////////// Attributes ///////////////////////////////////
	String myID;

	Coordinate home, work; // work/school or whatever

	public Status status = null;

	public int timeSinceDeparted = 0;

	public static int numActive = 0;
	String homeTract = "";
	String goalTract = "";
	String agentName = "";
	Node headquartersNode = null;
	Node goalNode = null;
	// private Point location;

	// private double basemoveRate = 20.0; // How much to move the agent by in each
	// step
	// private double moveRate = basemoveRate; // private double moveRate = 70;
	public Network familiarRoadNetwork = null;
	ArrayList<ArrayList<Edge>> familiarPaths = new ArrayList<ArrayList<Edge>>();

	Coordinate targetDestination = null;

	private double moveRate = 20.0;
	int minMoveRate = 20;
	int maxCarMoveRate = 60;
	int maxRoverMoveRate = 50;
	public double speed = 0; // useful for graph
	// private Color inboundColor = Color.black;
	// private Color outboundColor = Color.red;
	String type = "";

	double lastMove = -1;

	Coordinate goalCoord = null;

	////////// Movement + Goals //////////////////////////////
	double startIndex = 0.0; // start position of current line
	double endIndex = 0.0; // end position of current line
	double currentIndex = 0.0; // current location along line
	GeomPlanarGraphEdge currentEdge = null;
	ArrayList<Edge> pathFromHQToLSOA = new ArrayList<Edge>();
	int indexOnPath = 0;
	int pathDirection = 1;
	private LengthIndexedLine segment = null;
	PointMoveTo pointMoveTo = new PointMoveTo();

	public boolean getReachedGoal() {
		return hasResources;
	}

	public void setReachedGoal(boolean val) {
		hasResources = val;
	}

	int linkDirection = 1;

	////////// Resources ///////////////////////////////////
	public boolean hasResources = false;
	public static boolean active = false;
	public boolean activated = false;
	public int resources_Available = 0;
	public int resources_Distributed = 0;
	public int resources_Capacity = 0;

	////////// Parameters ///////////////////////////////////
	static private GeometryFactory fact = new GeometryFactory();
	public ArrayList<String> recordOfTrips = new ArrayList<String>();

	////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// AGENT /////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	public GeoNode getNode() {
		return node;
	}

	/**
	 * //////////////////////// Default Wrapper Constructor ///////////////////////
	 * Constructor: specifies parameters for the Agent's Default Wrapper Constructor
	 * 
	 * @param id
	 *            - unique string identifying the Agent
	 * @param position
	 *            - Coordinate indicating the initial position of the Agent
	 * @param home
	 *            - Coordinate indicating the Agent's home location
	 * @param work
	 *            - Coordinate indicating the Agent's workplace
	 * @param world
	 *            - reference to the containing Hotspots instance
	 */
	// public Agent(MK_7_1 g, String agentName, Coordinate homeTract, Coordinate
	// startingEdge, Coordinate goalEdge){
	// this(g, agentName, homeTract, startingEdge, goalEdge);
	// }

	/**
	 * //////////////////////// Specialised Constructor ////////////////////////////
	 * Constructor: used to specify parameters for an Agent
	 *
	 * @param location
	 *            - Coordinate indicating the initial position of the Agent
	 * @param headquartersNode
	 *            - Coordinate indicating the Agent's home location
	 * @param lsaoNode
	 *            - Coordinate indicating the Agent's LSOA destination
	 * @param world
	 *            - reference to the containing GloucestershireRouting instance
	 * 
	 *            Also Speed and various other parameters need adding
	 */
	// public MainAgent(MK_7_1 g, String agentName, Coordinate homeTract, Coordinate
	// startingEdge, Coordinate goalEdge) {
	public Agent(String id, Coordinate position, String homeTract, GeomPlanarGraphEdge startingEdge,
			GeomPlanarGraphEdge goalEdge) {

		super((new GeometryFactory()).createPoint(position)); // From Hotspots

		this.world = world;
		this.isMovable = true;
		Random random = null;

		// set up information about where the node is and where the Agent is going
		headquartersNode = startingEdge.getDirEdge(0).getFromNode();
		goalNode = goalEdge.getDirEdge(0).getToNode();
		this.agentName = agentName;
		this.homeTract = homeTract;
		this.goalTract = goalTract;

		// set the location to be displayed
		GeometryFactory fact = new GeometryFactory();

		// Now set up attributes for this agent
		if (this.world.random.nextBoolean()) {
			type = "4x4";
			this.addStringAttribute("TYPE", "4x4");
			System.out.println("Agent's Vehicle = " + type);

			int range = (int) (40.0 * this.world.random.nextGaussian());
			this.addIntegerAttribute("RANGE", range);
			// System. out.println("Agent's RANGE = " + range );

			resources_Available = (int) (Math.random() * 70) + 1;
			this.addIntegerAttribute("RESOURCES AVAILABLE", resources_Available);
			System.out.println("Agent's available resources = " + resources_Available);

			moveRate = Math.random() * ((maxRoverMoveRate - minMoveRate) + 1) + minMoveRate;
			// moveRate = (int) (Math.random() * 50) + 1;
			System.out.println("Agent's speed = " + moveRate);
			this.addDoubleAttribute("MOVE RATE", moveRate);
		} else {
			type = "Car";
			this.addStringAttribute("TYPE", "Car");
			System.out.println("Agent's Vehicle = " + type);

			int range = (int) (20.0 * this.world.random.nextGaussian());
			this.addIntegerAttribute("RANGE", range);
			// System. out.println("Agent's RANGE = " + range );

			resources_Available = (int) (Math.random() * 70) + 1;
			this.addIntegerAttribute("RESOURCES AVAILABLE", resources_Available);
			System.out.println("Agent's  available resources = " + resources_Available);

			moveRate = Math.random() * ((maxCarMoveRate - minMoveRate) + 1) + minMoveRate;
			// moveRate = (int) (Math.random() * 70) + 1;
			System.out.println("Agent's speed = " + moveRate);
			this.addDoubleAttribute("MOVE RATE", moveRate);
		}

		// Not everyone moves at the same speed
		// Assigns random speed
		// moveRate *= Math.abs(g.random.nextGaussian());
		// Assigns random speed between 0-70

		Coordinate startCoord = null;
		startCoord = headquartersNode.getCoordinate();
		Coordinate goalCoord = null;
		goalCoord = goalNode.getCoordinate();
		updatePosition(startCoord);
		// System.out.println("Agent's starting coordinate: " + this);
		System.out.println("Agent's starting coordinate: " + startCoord);
		System.out.println("Agent's goal coordinate: " + goalCoord);

		//////////////////////////////////////////////////////////////////////////////
		//////////////////////////// AGENT ATTRIBUTES ////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////

		/**
		 * ////////////////////////// Agent Resources ////////////////////////////////
		 *
		 * @return integer indicating the Agent's levels of resources
		 * 
		 */
		/*
		 * public resources(Integer resources_Available, Integer resources_Distributed,
		 * Integer resources_Capacity) { this.resources_Available = resources_Available;
		 * this.resources_Distributed = resources_Distributed; this.resources_Capacity =
		 * resources_Capacity; }
		 */

		/**
		 * ////////////////////////// Agent Colour ////////////////////////////////////
		 * Want to change the colour of the Agent's depending on their status: "heading
		 * back to HQ" or "heading to goal"
		 * 
		 */
		// public draw(Object object, Graphics2D graphics, DrawInfo2D info) { if(
		// hasResources ) graphics.setColor( outboundColor ); else graphics.setColor(
		// inboundColor );
		//
		// this code was copied from OvalPortrayal2D int x = (int)(info.draw.x -
		// info.draw.width / 20.0); int y = (int)(info.draw.y - info.draw.height /
		// 20.0); int width = (int)(info.draw.width); int height =
		// (int)(info.draw.height); graphics.fillOval(x,y,width, height); }
	}

	/**
	 * Based on the Agent's current activity and time of day, pick its next activity
	 */
	void pickDefaultActivity() {
		// if the Agent is moving, keep moving
		if (currentActivity == activity_travel && segment != null) {
			return;
		}
		// if the Agent is travelling but has no route,
		else if (currentActivity == activity_travel && segment == null) {

		}

	}

	/**
	 * Return the timestep that will correspond with the next instance of the given
	 * hour:minute combination
	 * 
	 * @param desiredHour
	 *            - the hour to find
	 * @param desiredMinuteBlock
	 *            - the minute to find
	 * @return the timestep of the next hour:minute combination
	 */
	int getTime(int desiredHour, int desiredMinuteBlock) {

		int result = 0;

		// the current time in the day
		int time = (int) (world.schedule.getTime());
		int numDaysSoFar = (int) Math.floor(time / 288);
		int currentTime = time % 288;

		int goalTime = desiredHour * 12 + desiredMinuteBlock;

		if (goalTime < currentTime)
			result = 288 * (numDaysSoFar + 1) + goalTime;
		else
			result = 288 * numDaysSoFar + goalTime;

		return result;
	}

	/**
	 * ////////////////////////// Step to Move Agent //////////////////////////////
	 * Called every tick by the scheduler. Moves the Agent along the path.
	 * 
	 * @param <osviInGeo>
	 */
	public void step(SimState state) {

		////////// Initial Checks //////////////////////////////

		// make sure the Agent is only being called once per tick

		if (lastMove >= state.schedule.getTime()) {
			System.out.println(this + "is being called too often!. " + "WHY??? EVERYBODY PANIC!!!");
			return;
		}

		// check that Agent has been placed on an Edge
		if (segment == null) {
			System.out.println(
					this + "'s segment == null. " + "Agent is NOT on an edge!!! " + "WHY??? EVERYBODY PANIC!!!");
			return;
		}

		// if(this.geometry.getCoordinate() == null) { System.out.println(this +
		// "'s location == null! " + "Agent has NO coordinate!!! " +
		// "WHY??? EVERYBODY PANIC!!!"); return; }

		// if (location.geometry.getCoordinate() == headquartersNode.getCoordinate()) {
		// System.out.println(location + " is at HQ ");
		// }

		if (this.geometry.getCoordinate() != goalNode.getCoordinate()
				|| this.geometry.getCoordinate() != headquartersNode.getCoordinate()) {
			// System.out.println(location + " is travelling.");
			distributing = false;
			inbound = false;
			status = Status.INBOUND;
			// System.out.println(location + " is " + status);
		}

		// check that we haven't already reached our destination
		else if (this.geometry.getCoordinate() == goalNode.getCoordinate()) {
			System.out.println(this + " IS AT ITS GOAL! ");
			status = Status.DISTRIBUTING;
			// setActive(gstate);
			System.out.println(this + " is " + status + this);

			//////////////////////////////////////////////////
			////// TODO CHANGE STATUS TO DISTRIBUTING //////
			////// DROP GOODS, DECREASE UNIT SCORE //////
			//////////////////////////////////////////////////

			recordOfTrips.add(
					this.toString() + " COMPLETED TRIP TO " + this.getGeometry().getCentroid().toString());
					//this.toString() + " COMPLETED TRIP TO " + this.getGeometry().geometry.getCentroid().toString());

			////////////////////////////////////////////////////////////////
			////// TODO RECORD VISITED GOAL AND MAKE IT AVAILABLE TO ///////
			////// getLargestUnassignedWard() ////////
			////////////////////////////////////////////////////////////////

			flipPath();
			state.schedule.scheduleOnce(state.schedule.getTime() + 50, this);
			// makes Agent wait before flipping route

			return;
		}

		// make sure that we're heading in the right direction
		// boolean toWork = ((MK_7) state).goToWork;
		// if ((toWork && pathDirection < 0) || (!toWork && pathDirection > 0)) {
		// flipPath();
		// }
		speed = progress(moveRate);
		currentIndex += speed;

		// check to see if the progress has taken the current index beyond its goal
		// given the direction of movement. If so, proceed to the next edge
		if (linkDirection == 1 && currentIndex > endIndex) {
			Coordinate currentPos = segment.extractPoint(endIndex);
			updatePosition(currentPos);
			transitionToNextEdge(currentIndex - endIndex);
		} else if (linkDirection == -1 && currentIndex < startIndex) {
			Coordinate currentPos = segment.extractPoint(startIndex);
			updatePosition(currentPos);
			transitionToNextEdge(startIndex - currentIndex);
		} else {
			// just update the position!
			Coordinate currentPos = segment.extractPoint(currentIndex);

			updatePosition(currentPos);
			// System.out.println(location + " currentPos " + currentPos);
		}

		state.schedule.scheduleOnce(this);
		// lastMove = state.schedule.getTime();
	}

	////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// METHODS ///////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	/**
	 * ////////////////////////// Set Up Path /////////////////////////////////////
	 * a course to take the Agent to the given coordinates
	 *
	 * @param place
	 *            - the target destination
	 * @return 1 for success, -1 for a failure to find a path, -2 for failure based
	 *         on the provided destination or current position
	 */
	// int headFor(Coordinate place, Network roadNetwork) {

	/*
	 * first, record from where the agent is starting
	 * 
	 * if the current node and the current edge don't match, there's a problem with
	 * the Agent's understanding of its current position
	 * 
	 * ///////////////////// FINDING THE GOAL //////////////////
	 * 
	 * set up goal information
	 * 
	 * be sure that if the target location is not a node but rather a point along an
	 * edge, that point is recorded
	 * 
	 * ///////////////////// FINDING A PATH /////////////////////
	 * 
	 * if it fails, give up
	 * 
	 * /////////////////////CHECK FOR BEGINNING OF PATH ////////
	 * 
	 * we want to be sure that we're situated on the path *right now*, and that if
	 * the path doesn't include the link we're on at this moment that we're both a)
	 * on a link that connects to the startNode b) pointed toward that startNode
	 * Then, we want to clean up by getting rid of the edge on which we're already
	 * located Make sure we're in the right place, and face the right direction
	 * 
	 * reset stuff
	 * 
	 * /////////////////////CHECK FOR END OF PATH //////////////
	 * 
	 * we want to be sure that if the goal point exists and the Agent isn't already
	 * on the edge that contains it, the edge that it's on is included in the path
	 * 
	 * make sure the point is on the last edge
	 * 
	 * set up the coordinates
	 */

	////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// end METHODS ///////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// ROUTING ///////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////

	/**
	 * ////////////////////////// A* Route Initialisation /////////////////////////
	 * Initialization of an Agent: find an A* path from HQ to LSOA!
	 *
	 * @param state
	 * @return whether or not the Agent successfully found a path to work
	 */
	/*public boolean start(MK_7_1 state) {
		findNewAStarPath(state);
		if (pathFromHQToLSOA.isEmpty()) {
			System.out.println("Initialization of Agent (" + homeTract
					+ ") failed: it is located in a part of the network that cannot " + "access the given goal.");
			return false;
		} else {
			System.out.println("Agent has been assigned a path.");
			return true;
		}
	}
	 */
	/**
	 * Set up a course to take the Agent to the given coordinates
	 * 
	 * @param place
	 *            - the target destination
	 * @return 1 for success, -1 for a failure to find a path, -2 for failure based
	 *         on the provided destination or current position
	 */
	int headFor(Coordinate place, Network roadNetwork) {

		// first, record from where the agent is starting
		startPoint = this.geometry.getCoordinate();
		goalPoint = null;

		// if the current node and the current edge don't match, there's a problem with
		// the Agent's understanding of its
		// current position
		if (!(edge.getTo().equals(node) || edge.getFrom().equals(node))) {
			System.out.println((int) world.schedule.getTime() + "\t" + this.myID
					+ "\tMOVE_ERROR_mismatch_between_current_edge_and_node");
			return -2;
		}

		// FINDING THE GOAL //////////////////

		// set up goal information
		targetDestination = this.snapPointToRoadNetwork(place);

		GeoNode destinationNode = world.getClosestGeoNode(targetDestination);// place);
		if (destinationNode == null) {
			System.out.println(
					(int) world.schedule.getTime() + "\t" + this.myID + "\tMOVE_ERROR_invalid_destination_node");
			return -2;
		}

		// be sure that if the target location is not a node but rather a point along an
		// edge, that
		// point is recorded
		if (destinationNode.geometry.getCoordinate().distance(targetDestination) > MK_7_1.resolution)
			goalPoint = targetDestination;
		else
			goalPoint = null;

		// FINDING A PATH /////////////////////

		if (path == null)
			path = pathfinder.astarPath(node, destinationNode, roadNetwork);

		// if it fails, give up
		if (path == null) {
			return -1;
		}

		// CHECK FOR BEGINNING OF PATH ////////

		// we want to be sure that we're situated on the path *right now*, and that if
		// the path
		// doesn't include the link we're on at this moment that we're both
		// a) on a link that connects to the startNode
		// b) pointed toward that startNode
		// Then, we want to clean up by getting rid of the edge on which we're already
		// located

		// Make sure we're in the right place, and face the right direction
		if (edge.getTo().equals(node))
			direction = 1;
		else if (edge.getFrom().equals(node))
			direction = -1;
		else {
			System.out.println((int) world.schedule.getTime() + "\t" + this.myID
					+ "MOVE_ERROR_mismatch_between_current_edge_and_node_2");
			return -2;
		}

		// reset stuff
		if (path.size() == 0 && targetDestination.distance(geometry.getCoordinate()) > world.resolution) {
			path.add(edge);
			node = (GeoNode) edge.getOtherNode(node); // because it will look for the other side in the navigation!!!
														// Tricky!!
		}

		// CHECK FOR END OF PATH //////////////

		// we want to be sure that if the goal point exists and the Agent isn't already
		// on the edge
		// that contains it, the edge that it's on is included in the path
		if (goalPoint != null) {// && path.size() > 0) {

			ListEdge myLastEdge = world.getClosestEdge(goalPoint);

			if (myLastEdge == null) {
				System.out.println((int) world.schedule.getTime() + "\t" + this.myID
						+ "\tMOVE_ERROR_goal_point_is_too_far_from_any_edge");
				return -2;
			}

			// make sure the point is on the last edge
			Edge lastEdge;
			if (path.size() > 0)
				lastEdge = path.get(0);
			else
				lastEdge = edge;

			Point goalPointGeometry = world.fa.createPoint(goalPoint);
			if (!lastEdge.equals(myLastEdge)
					&& ((MasonGeometry) lastEdge.info).geometry.distance(goalPointGeometry) > MK_7_1.resolution) {
				if (lastEdge.getFrom().equals(myLastEdge.getFrom()) || lastEdge.getFrom().equals(myLastEdge.getTo())
						|| lastEdge.getTo().equals(myLastEdge.getFrom()) || lastEdge.getTo().equals(myLastEdge.getTo()))
					path.add(0, myLastEdge);
				else {
					System.out.println((int) world.schedule.getTime() + "\t" + this.myID
							+ "\tMOVE_ERROR_goal_point_edge_is_not_included_in_the_path");
					return -2;
				}
			}

		}

		// set up the coordinates
		this.startIndex = segment.getStartIndex();
		this.endIndex = segment.getEndIndex();

		return 1;
	}

	/**
	 * ////////////////////////// Plot A* Route ///////////////////////////////////
	 * Plots a path between the Agent's home Node (HQ) and its goal Node (LSOA)
	 */
	
	private void findNewAStarPath(MK_7_1 geoTest) {

		// get the home and work Nodes with which this Agent is associated
		Node currentJunction = geoTest.roadNetwork.findNode(this.geometry.getCoordinate());
		Node destinationJunction = goalNode;

		if (currentJunction == null) {
			return; // just a check
		}
		// find the appropriate A* path between them
		AStar pathfinder = new AStar();
		ArrayList<Edge> path = pathfinder.astarPath(currentJunction, destinationJunction, geoTest.roadNetwork);

		// if the path works, lay it in
		if (path != null && path.size() > 0) {

			// save it
			pathFromHQToLSOA = path;

			// set up how to traverse this first link
			GeomPlanarGraphEdge edge = (GeomPlanarGraphEdge) path.get(0).getEdge();
			setupEdge(edge);

			// update the current position for this link
			updatePosition(segment.extractPoint(currentIndex));
		}
	}

	double progress(double val) {
		double edgeLength = ((GeomPlanarGraphEdge) currentEdge).getLine().getLength();
		double traffic = world.edgeTraffic.get(currentEdge).size();
		double factor = 1000 * edgeLength / (traffic * 5);
		factor = Math.min(1, factor);
		return val * linkDirection * factor;
	}
 	
	/**
	 * ////////////////////// Number of Active Agents ////////////////////////////
	 * 
	 */
	private void setActiveNumAgents(MK_7_1 gState) {
		boolean alreadyActive = false;

		if (Agent.active)
			alreadyActive = true;

		int numActive = 0;

		// Loop through all the agents within vision and count the number of active
		// civilians
		/*
		 * MasonGeometry buffer = new
		 * MasonGeometry(this.location.buffer(osviState.personVision), this); Bag
		 * persons = osviState.persons.getCoveredObjects(buffer); Bag cops =
		 * osviState.cops.getCoveredObjects(buffer); for(Object person : MainAgent){
		 * MainAgent p = (MainAgent)((MasonGeometry) person).getUserData();
		 */
		if (agents.Agent.active) {
			numActive++;
		}
		// }

		/* Calculations for going active */
		// arrestProbability = 1 - Math.exp(-2.3*((cops.size()/(numActive+1))));
		// grievance = perceivedHardship * (1 - govtLegitimacy);
		// this.active = (grievance - (riskAversion * arrestProbability)) >
		// osviState.threshold;

		if (!alreadyActive && this.active) {
			gState.activeCount++;
			activated = true;
		} else if (alreadyActive && !this.active) {
			gState.activeCount--;
		}
	}

	/**
	 * ////////////////////////// Flip Agent's Route /////////////////////////////
	 * Flip the Agent's path around
	 */
	public void flipPath() {
		distributing = false;
		inbound = false;
		status = Status.INBOUND;
		System.out.println(this + " is " + status);
		// this.timeSinceDeparted = 0;

		pathDirection = -pathDirection;
		linkDirection = -linkDirection;
	}

	/**
	 * ////////////////////////// Move Agent to Next Edge ////////////////////////
	 * Transition to the next edge in the path
	 * 
	 * @param residualMove
	 *            - the amount of distance the Agent can still travel this turn
	 */
	
	void transitionToNextEdge(double residualMove) {

		// update the counter for where the index on the path is
		indexOnPath += pathDirection; // indexOnPath = indexOnPath + pathDirection

		// check to make sure the Agent has not reached the end
		// of the path already
		if ((pathDirection > 0 && indexOnPath >= pathFromHQToLSOA.size()) || (pathDirection < 0 && indexOnPath < 0))
		// if pathDirection GREATER 0 AND indexOnPath EQUAL or GREATER
		// pathFromHQToLSOA.size OR pathDirection LESS 0 AND indexOnPath LESS 0
		// depends on where you're going!
		{
			status = Status.OUTBOUND;
			distributing = true;
			inbound = true;
			// System.out.println(this + " is " + status);
			indexOnPath -= pathDirection; // make sure index is correct
			return;
		}

		// move to the next edge in the path
		GeomPlanarGraphEdge edge = (GeomPlanarGraphEdge) pathFromHQToLSOA.get(indexOnPath).getEdge();
		setupEdge(edge);
		speed = progress(residualMove);
		currentIndex += speed;

		// check to see if the progress has taken the current index beyond its goal
		// given the direction of movement. If so, proceed to the next edge
		if (linkDirection == 1 && currentIndex > endIndex) {
			transitionToNextEdge(currentIndex - endIndex);
		} else if (linkDirection == -1 && currentIndex < startIndex) {
			transitionToNextEdge(startIndex - currentIndex);
		}
	}
	 
	
	/**
	 * ////////////////////////// Agent's Route Info /////////////////////////////
	 * Sets the Agent up to proceed along an Edge
	 * 
	 * @param edge
	 *            - the GeomPlanarGraphEdge to traverse next
	 */
	void setupEdge(GeomPlanarGraphEdge edge) {

		// clean up on old edge
		if (currentEdge != null) {
			ArrayList<Agent> traffic = world.edgeTraffic.get(currentEdge);
			traffic.remove(this);
		}
		currentEdge = edge;

		// update new edge traffic
		if (world.edgeTraffic.get(currentEdge) == null) {
			world.edgeTraffic.put(currentEdge, new ArrayList<Agent>());
		}
		world.edgeTraffic.get(currentEdge).add(this);

		// set up the new segment and index info for the Agent's position on the road
		// segment
		LineString line = ((GeomPlanarGraphEdge) edge).getLine();
		segment = new LengthIndexedLine(line);
		// segment = new
		// LengthIndexedLine((LineString)((MasonGeometry)edge.info).geometry); // From
		// Hotspots
		startIndex = segment.getStartIndex();
		endIndex = segment.getEndIndex();
		// currentIndex = segment.indexOf(position); // From Hotspots
		linkDirection = 1;

		// check to ensure that Agent is moving in the right direction
		double distanceToStart = line.getStartPoint().distance(this.geometry),
				distanceToEnd = line.getEndPoint().distance(this.geometry);
		if (distanceToStart <= distanceToEnd) { // closer to start
			currentIndex = startIndex;
			linkDirection = 1;
		} else if (distanceToEnd < distanceToStart) { // closer to end
			currentIndex = endIndex;
			linkDirection = -1;
		}
	}

	/**
	 * Snap the coordinate to the nearest point on the road network
	 * 
	 * @param c
	 *            - the point in question
	 * @return - nearest point on the road network
	 */
	public Coordinate snapPointToRoadNetwork(Coordinate c) {
		network.ListEdge myEdge = null;
		double resolution = MK_7_1.resolution;

		// if the network hasn't been properly set up, don't try to find something on it
		// =\
		if (world.networkEdgeLayer.getGeometries().size() == 0)
			return null;

		// while there's no edge, expand outward until the Agent finds one
		while (myEdge == null) {
			myEdge = world.getClosestEdge(c, resolution);
			resolution *= 10;
		}

		// having found a line, find the index of the point on that line
		LengthIndexedLine closestLine = new LengthIndexedLine(
				(LineString) (((MasonGeometry) myEdge.info).getGeometry()));
		double myIndex = closestLine.indexOf(c);
		return closestLine.extractPoint(myIndex);
	}

	/**
	 * ////////////////////////// Move Agent /////////////////////////////////////
	 * Move the Agent to the given coordinates
	 */
	public void updatePosition(Coordinate c) {
		PointMoveTo p = new PointMoveTo();
		p.setCoordinate(c);
		geometry.apply(p);
		geometry.geometryChanged();
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// end ROUTING /////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// UTILITIES ///////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	// ///////////////////// GETTERS /////////////////////
	public Coordinate getHome() {
		return home;
	}

	public Coordinate getWork() {
		return work;
	}

	public int getActivity() {
		return this.currentActivity;
	}
	// public double getValence(){ return this.stress; }

	/**
	 * ////////////////////////// Agent's Location //////////////////////////////
	 * Return geometry representing Agent's location
	 */
	// public MasonGeometry getGeometry() {
	// return location;
	// }

	/**
	 * ////////////////////////// Default Path //////////////////////////////////
	 * Generates default paths between home and work, to ease on computation during
	 * the actual run
	 */
	public void setupPaths() {
		if (work != null) {
			GeoNode workNode = world.getClosestGeoNode(this.work);
			GeoNode homeNode = world.getClosestGeoNode(this.home);

			ArrayList<Edge> pathFromHomeToWork = pathfinder.astarPath(homeNode, workNode, world.roadNetwork);
			this.familiarPaths.add(pathFromHomeToWork);

			ArrayList<Edge> pathFromWorkToHome = pathfinder.astarPath(workNode, homeNode, world.roadNetwork);
			this.familiarPaths.add(pathFromWorkToHome);
		}
	}

	/**
	 * ////////////////////////// Step Wrapper //////////////////////////////////
	 * Wrapper around step, so that it can be called from other functions
	 */
	void stepWrapper() {
		this.step(world);
	}

	public boolean start(MK_7_1 mk_7_1) {
		// TODO Auto-generated method stub
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////// end UTILITIES ///////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

}