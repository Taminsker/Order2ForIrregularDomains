#include "funtovec.h"

FunToVec::FunToVec (Mesh * mesh) :
    m_mesh (mesh)
{}

FunToVec::~FunToVec ()
{}

Vector FunToVec::Make (double (*f)(Point))
{
    /// - Créer un vecteur.
    /// - Boucler sur les points du maillage.
    /// - Appliquer la fonction f à chaque points.
}

Vector FunToVec::Make (double (*f)(Point, double), double t)
{
    /// - Créer un vecteur.
    /// - Boucler sur les points du maillage.
    /// - Appliquer la fonction f à chaque points et au temps t.
}

Vector FunToVec::Make (double value)
{
    /// - Créer un vecteur.
    /// - Initialiser le vecteur à value.
}
