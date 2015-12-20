import java.util.ArrayList;

public class Edge
{
	public Vertex a;
	public Vertex b;
	public int direction;
	public ArrayList<Vertex> edgeSquares;
	
	public Edge(Vertex a, Vertex b, int direction)
	{
		this.a = a;
		this.b = b;
		this.direction = direction;
		edgeSquares = new ArrayList<Vertex>();
		addLine(a.y, a.x, b.y, b.x);
	}
	
	public void addSquare(Vertex v)
	{
		edgeSquares.add(v);
	}
	
	public void removeSquare(Vertex v)
	{
		edgeSquares.remove(v);
	}
	
	public void addLine(int y1, int x1, int y2, int x2)
	{
		if (y1 == y2)
		{
			while (x1 > x2)
			{
				x1--;
				addSquare(new Vertex(y1, x1, 0));
			}
			while (x1 < x2)
			{
				x1++;
				addSquare(new Vertex(y1, x1, 0));
			}
		}
		else
		{
			while (y1 > y2)
			{
				y1--;
				addSquare(new Vertex(y1, x1, 0));
			}
			while (y1 < y2)
			{
				y1++;
				addSquare(new Vertex(y1, x1, 0));
			}
		}			
	}
	
	public String toString()
	{
		return "E: " + a.x + "," + a.y + " - " + b.x + "," + b.y;
	}
}