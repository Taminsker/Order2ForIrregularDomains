#include "writer.h"

Writer::Writer (Mesh * mesh) :
    m_mesh (mesh)
{

}

Writer::~Writer ()
{

}
void Writer::SetVectorAnalytical (Vector *vector)
{
    m_sol_ana = vector;
    return;
}

void Writer::SetVectorNumerical (Vector *vector)
{
    m_sol_num = vector;
    return;
}

void Writer::SetVectorErrorAbs (Vector *vector)
{
    m_error_abs = vector;
    return;
}

void Writer::SetFilename (std::string filename)
{
    m_filename = filename;
    return;
}

void Writer::WriteBothDomainsOn ()
{
    m_domain = true;
    return;
}

void Writer::WriteBothDomainsOff ()
{
    m_domain = false;
    return;
}

void Writer::SetCurrentIteration (int i)
{
    m_index = i;
    return;
}

void Writer::WriteNow ()
{
    return;
}
