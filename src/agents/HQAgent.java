package agents;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.ArrayList;

import network.AStar;
import sim.MK_7_1;
import sim.Status;
import sim.app.geo.sickStudents.SickStudentsModel;
import sim.app.geo.sickStudents.Student;
import sim.app.geo.sickStudents.School.SchoolType;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.portrayal.DrawInfo2D;
import sim.util.geo.GeomPlanarGraphDirectedEdge;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import sim.util.geo.PointMoveTo;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.linearref.LengthIndexedLine;
import com.vividsolutions.jts.planargraph.Node;

/**
*
* "MK_7_" is iteration 7.1 of my EngD project model. It varies little from MK_7, 
* and just has updated details to include Gloucestershire date. It was originally
* adapted from the MASON demo, "Gridlock", made by Sarah Wise, Mark Coletti, and 
* Andrew Crooks.
* 
* The model reads a number of GIS shapefiles and displays a road network, two 
* Environment Agency flood maps and a bespoke Open Source Vulnerability Index 
* (OSVI). The model reads in a .CSV and generates a predetermined number of agents 
* with set characteristics. The agents are placed on the road network and are 
* located at a Red Cross office. The model reads a separate .CSV and assigns goal 
* locations to each agent at random from a predetermined list. The agents are 
* assigned speeds at random. Once the model is started, the agents move from 
* A to B, then they change direction and head back to their start position. 
* The process repeats until the user quits.
*
* @author KJGarbutt
*
*/
@SuppressWarnings("serial")
public class HQAgent implements Steppable	{
    
		MK_7_1 world;
    	ArrayList<agents.MainAgent> agentList = new ArrayList<agents.MainAgent>();
    	public String name;
    	public SchoolType type;
    	public boolean closed = false;
    	
    	public int catchmentCount = 0;
    	
    	private Student getRandomStudent(Student butNotThisStudent) {
    		Student s;
    		do
            {
                s = MainAgent.get(world.random.nextInt(MainAgent.size()));
            } while (!(!s.homebound && (s != butNotThisStudent)));
    		
    		return s;
    	}
    	
    	public double getProportionOfSickStudents() {
    		if (MainAgent.isEmpty())
            {
                return 0;
            }
    		
    		int sick = 0;
    		for (Student s : MainAgent)
            {
                if (s.status == Status.INFECTED)
                {
                    sick++;
                }
            }
    		
    		return sick / (double)MainAgent.size();
    	}
    	
    	public double getProportionOfInboundAgents() {
    		if (MainAgent.isEmpty())
            {
                return 0;
            }
    		
    		int homebound = 0;
    		for (Student s : MainAgent)
            {
                if (s.homebound)
                {
                    homebound++;
                }
            }
    		
    		return homebound / (double)MainAgent.size();
    	}
    	
    	@Override
    	public void step(SimState state) {
    		if (closed)
            {
                return;
            }
    		
    		int inAttendence = 0;
    		for (Student s : MainAgent)
            {
                if (!s.homebound)
                {
                    inAttendence++;
                }
            }
    		
    		if (inAttendence < 2)
            {
                return;
            }
    	}
    }
