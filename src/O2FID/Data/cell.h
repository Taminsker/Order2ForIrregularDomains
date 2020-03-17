#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <vector>

typedef enum {
    IN_DOMAIN_INTERN_OMEGA,
    IN_DOMAIN_EXTERN_OMEGA,
    IN_MIX_DOMAIN
} CELL_LOCATION;

class Point;

class Cell
{
public:
    Cell ();
    ~Cell ();

    void AddPoint (Point * p);
    void RemovePoint (Point * p);
    CELL_LOCATION GetLocate () const;

    int GetType () const;

    int GetNumberOfInfos () const;

    friend std::ostream & operator<< (std::ostream &out, const Cell &c);

protected:
    std::vector <Point *> m_points;

};
std::ostream & operator<< (std::ostream &out, const Cell &c);


#endif // CELL_H
