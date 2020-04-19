#include "funtovec.h"

int NumberOfCall = 0;

Vector FunToVec (Mesh * mesh, double (*f) (Point, double), double t)
{
    /* Nombre de point total dans le mesh à calculer. */
    int n = mesh->GetNumberOfTotalPoints ();
    /* Création d'un vecteur de n points. */
    Vector Fvec (n);

    for (int i = 0; i < n; i++)
    {
        /* On récupère le i-ème point du maillages. */
        Point p_i = *mesh->GetPoint (i);

        /* On affecte la valeur de f à ce point. */
        Fvec (i) = f(p_i, t);
    }

    return Fvec;
}

Vector FunToVec (Mesh * mesh, double value)
{
    //    /* Nombre de point total dans le mesh à calculer. */
    int n = mesh->GetNumberOfTotalPoints ();
    //    /* Création d'un vecteur de n points. */
    Vector Fvec (n);

    //    for (int i = 0; i < n; i++)
    //    {
    //        Fvec (i) = value;
    //    }

    return Fvec.setConstant (value);
}


Vector TestVec (Mesh * mesh)
{

    int numPoints = mesh->GetNumberOfTotalPoints ();
    // Vecteur de sortie
    Vector out (numPoints);
    out.setConstant (0);

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    int compt = 0;
    int i = 0;
    int j = 0;
    int k = 0;

    std::cout << "Nx " << Nx << std::endl;
    std::cout << "Ny " << Ny << std::endl;
    std::cout << "Nz " << Nz << std::endl;

    int N = int(std::max(Nx, std::max(Ny, Nz)) / 2) +1;
    std::cout << "N " << N << std::endl;

    //    for (int d = 0; d < N; ++d)
    //    {

    //        j = d;
    //        for (i = d; i < Nx-d; ++i)
    //        {
    //            int index = (k * Ny + j) * Nx + i;
    //            printf ("(1) index : %i/%i (%i, %i, %i) --> %i\n", index, numPoints, i, j, k, compt);
    //            out.coeffRef (index) = compt;
    //            compt++;
    //        }

    //        i = Nx-1 - d;
    //        for (j = d+1; j < Ny-1-d; ++j)
    //        {
    //            int index = (k * Ny + j) * Nx + i;
    //            printf ("(2) index : %i/%i (%i, %i, %i) --> %i\n", index, numPoints, i, j, k, compt);
    //            out.coeffRef (index) = compt;
    //            compt++;
    //        }

    //        if(Ny-1-d != d)
    //        {
    //            j = Ny-1-d;
    //            for (i = Nx-1-d; i >= d; --i)
    //            {
    //                int index = (k * Ny + j) * Nx + i;
    //                printf ("(3) index : %i/%i (%i, %i, %i) --> %i\n", index, numPoints, i, j, k, compt);
    //                out.coeffRef (index) = compt;
    //                compt++;
    //            }
    //        }

    //        if (d != Nx-1-d)
    //        {
    //            i = d;
    //            for (j = Ny-1-1-d; j >= d+1; --j)
    //            {
    //                int index = (k * Ny + j) * Nx + i;
    //                printf ("(4) index : %i/%i (%i, %i, %i) --> %i\n", index, numPoints, i, j, k, compt);
    //                out.coeffRef (index) = compt;
    //                compt++;
    //            }
    //        }
    //    }

    for (int k = 0; k < Nz; ++k)
    {
        for (int d = 0; d < N; ++d)
        {
            printf("New loop on d : %i\n", d);

            std::vector<std::vector<std::vector<int>>> v;


            if (d < Nx-d)
                v.push_back (Build ({d, Nx-d-1}, d, k));

            if (d+1<Ny-1-d)
                v.push_back (Build (Nx-1-d, {d+1, Ny-d-2}, k));

            if(Ny-1-d != d && Nx-1-d >= d)
                v.push_back (Build ({Nx-1-d, d}, Ny-1-d, k));
            if (d != Nx-1-d && Ny-2-d >= d+1)
                v.push_back (Build (d, {Ny-2-d, d+1}, k));

            for (auto a : v)
            {
                printf("\tNew loop\n");

                auto i_vec = a.at (0);
                auto j_vec = a.at (1);
                auto k_vec = a.at (2);

                for (size_t c = 0; c < i_vec.size (); ++c)
                {
                    int i = i_vec.at (c);
                    int j = j_vec.at (c);
                    int k = k_vec.at (c);


                    int index = (k * Ny + j) * Nx + i;
                    printf ("\t\tIndex : %i/%i (%i, %i, %i) --> %i\n", index, numPoints, i, j, k, compt);
                    out.coeffRef (index) = compt;
                    compt++;
                }


            }
        }
    }

    return out;

}

Vector FunToVec (Mesh * mesh,
                 Point (*f) (double, double),
                 Point (*fprim) (double, double),
                 size_t numOfPtsBorder,
                 double t,
                 double a,
                 double b)
{
    int numPoints = mesh->GetNumberOfTotalPoints ();
    // Vecteur de sortie
    Vector out (numPoints);
    out.setConstant (-50);

    double h = double (b - a) / double (numOfPtsBorder);

    std::vector<Point *> position (numOfPtsBorder);

    for (size_t i = 0; i < numOfPtsBorder; ++i)
            position.at (i) = new Point (f (a + i * h, t));

    double meanDist = 0.;
    for (size_t i = 0; i < numOfPtsBorder - 1; ++i)
        meanDist += EuclidianDist (*(position.at (i)), *(position.at (i+1)));

    meanDist /= numOfPtsBorder - 1;

    double ratio = numOfPtsBorder * 5 * 1e-2 / meanDist;

    int middle = int (numOfPtsBorder / 2) + 1;
    int area = middle;
    int lastIdx = 0;

    Point * p_last = mesh->GetPoint (0);

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    int N = int(std::max(Nx, std::max(Ny, Nz)) / 2) +1;

    for (int k = 0; k < Nz; ++k)
    {
        for (int d = 0; d < N; ++d)
        {
//            printf("New loop on d : %i\n", d);

            std::vector<std::vector<std::vector<int>>> vec;


            if (d < Nx-d)
                vec.push_back (Build ({d, Nx-d-1}, d, k));

            if (d+1<Ny-1-d)
                vec.push_back (Build (Nx-1-d, {d+1, Ny-d-2}, k));

            if(Ny-1-d != d && Nx-1-d >= d)
                vec.push_back (Build ({Nx-1-d, d}, Ny-1-d, k));
            if (d != Nx-1-d && Ny-2-d >= d+1)
                vec.push_back (Build (d, {Ny-2-d, d+1}, k));

            for (auto empla : vec)
            {
//                printf("\tNew loop\n");

                auto i_vec = empla.at (0);
                auto j_vec = empla.at (1);

                for (size_t empla_in = 0; empla_in < i_vec.size (); ++empla_in)
                {
                    int i = i_vec.at (empla_in);
                    int j = j_vec.at (empla_in);


                    Point* p = mesh->GetPoint (i, j, k);

                    int idx =p->GetGlobalIndex ();

                    if (p_last != p)
                        area = int (EuclidianDist (*p, *p_last) * ratio);

//                    printf("\t area: %i, middle: %i\n", area, middle);
                    area = std::min (middle, area);

                    while (true)
                    {
                        //            printf ("MS - (%i)p : %i, lastIdx : %i, area : %i\n", idx, p->GetGlobalIndex (), lastIdx, area);
                        int idxMin = 0;
                        double minima_dist = numOfPtsBorder;

                        int first = Remainder (lastIdx - area, int (numOfPtsBorder));

                        int last = Remainder (lastIdx + area, int (numOfPtsBorder));

                        //            printf ("\tfirst : %i, last : %i\n", lastIdx - area, lastIdx + area);
                        //            printf ("\tfirst : %i, last : %i\n", first, last);
                        for (int i = -area; i <= area; ++i)
                        {
                            int idxCurve = Remainder (lastIdx - i, int(numOfPtsBorder));
                            //                printf ("\ti : %i, idxCurve : %i\n", i, idxCurve);

                            double dist = EuclidianDist (*p, f(a + idxCurve * h, t));

                            if (dist < minima_dist)
                            {
                                idxMin = idxCurve;
                                minima_dist = dist;
                            }

                            //                printf ("\ti : %i, idxCurve : %i idxMin : %i dist %f\n", i, idxCurve, idxMin, dist);

                        }

                        lastIdx = idxMin;

                        //            printf("\t\tlastIdx : %i, first %i, last %i\n", lastIdx, first, last);

                        if (lastIdx != first && lastIdx != last)
                            break;
                    }

//                    printf ("MS - (%i)p : %i, lastIdx : %i, area : %i\n", idx, p->GetGlobalIndex (), lastIdx, area);


                    Point p_tangent = fprim (a + lastIdx * h, t);
                    Point normal_out = {-p_tangent.y, p_tangent.x, 0.};

                    out.coeffRef (idx) = (normal_out | (*(position.at (size_t (lastIdx)))- *p));

                    p_last = p;
                }
            }
        }
    }

    // Suppression du vecteur de position
    for (Point * p : position)
        delete p;

    position.clear ();

    // Fin de l'algo
    return out;
}



//Vector FunToVec (Mesh * mesh,
//                 Point (*f) (double, double),
//                 Point (*fprim) (double, double),
//                 double a,
//                 double b,
//                 size_t numOfPtsBorder, double t)
//{

//    WagonF2V wgn;
//    wgn.size = numOfPtsBorder;


//    // Vecteur de sortie
//    int N = mesh->GetNumberOfTotalPoints ();
//    Vector outVector (N);
//    outVector.setConstant (-50);

//    wgn.out = &outVector;


//    // Vecteurs position et vitesse
//    double hOnBorder = double (b - a) / double (wgn.size);
//    std::vector<Point *> position (wgn.size);
//    std::vector<Point *> velocity (wgn.size);
//    std::vector<bool> reached (size_t (N), false);

//    for (size_t i = 0; i < wgn.size; ++i)
//    {
//        position.at (i) = new Point (f (a + i * hOnBorder, t));
//        velocity.at (i) = new Point (fprim (a + i * hOnBorder, t));

//        position.at (i)->SetGlobalIndex (int(i));
//        velocity.at (i)->SetGlobalIndex (int(i));
//    }

//    wgn.position = &position;
//    wgn.velocity = &velocity;
//    wgn.reached = &reached;

//    // DistanceMoyenne
//    double meanDist = 0.;
//    for (int i = 0; i < int(position.size () - 1); ++i)
//    {
//        auto a = position.at (size_t (i));
//        //        printf("\t\n b %i \n", a->GetGlobalIndex ());

//        meanDist += EuclidianDist (*(position.at (size_t (i))),
//                                   *(position.at(size_t (i+1))));
//    }


//    meanDist /= (position.size () - 1);

//    double ratio = numOfPtsBorder * 1.2 * 1e-1 / meanDist;

//    wgn.ratio = ratio;

//    // Appel de la fonction Spreading
//    SpreadingFrom (&wgn, mesh->GetPoint (0));

//    // Suppression du vecteur de position
//    for (Point * p : position)
//        delete p;

//    position.clear ();

//    // Suppression du vecteur vitesse
//    for (Point * p : velocity)
//        delete p;

//    velocity.clear ();

//    // Fin de l'algo
//    return outVector;

//}

//void SpreadingFrom (WagonF2V * wgn,
//                    Point * p,
//                    double lastDst,
//                    int lastIdx)
//{
//    int p_index = p->GetGlobalIndex ();

//    if (wgn->reached->at (size_t (p_index)))
//        return;

//    double* coeff = &(wgn->out->coeffRef (p_index));

//    if (lastIdx < 0)
//        lastIdx = MinimumSearch (wgn, coeff, p, 0, int(wgn->size));
//    else
//        lastIdx = MinimumSearch (wgn, coeff, p, lastIdx, int (wgn->ratio * lastDst));

//    wgn->reached->at (size_t (p_index)) = true;

//    for (Point * p_neigh : p->GetListNeighbours ())
//        SpreadingFrom (wgn, p_neigh, EuclidianDist(*p, *p_neigh), lastIdx);

//    return;
//}

//int MinimumSearch (WagonF2V* wgn,
//                   double * coeff,
//                   Point* p,
//                   int lastIdx,
//                   int area)
//{
//    NumberOfCall++;
//    printf("Call : %i\n", NumberOfCall);
//    area = std::min (int(wgn->size / 2) + 1, area);

//    //    printf ("Call MS for p : %i, area : %i, lastIndex : %i, max : %i\n", p->GetGlobalIndex (), area, lastIdx, wgn->size);

//    size_t idxOfMin = 0;
//    double minima_dist = wgn->size;

//    for (int i = -area; i <= area; ++i)
//    {
//        int index = lastIdx - i;

//        while (index >= int(wgn->size) || index < 0)
//            index += (index < 0 ? 1 : -1) * int(wgn->size);

//        Point* p_curve = wgn->position->at (size_t (index));

//        double dist = EuclidianDist (*p, *p_curve);

//        if (dist < minima_dist)
//        {
//            idxOfMin = size_t (index);
//            minima_dist = dist;
//        }
//    }

//    int first = lastIdx - area;
//    while (first >= int(wgn->size) || first < 0)
//        first += (first < 0 ? 1 : -1) * int(wgn->size);

//    int last = lastIdx + area;
//    while (last >= int(wgn->size) || last < 0)
//        last += (last < 0 ? 1 : -1) * int(wgn->size);

//    if (idxOfMin == size_t (first) || idxOfMin == size_t (last))
//        return MinimumSearch (wgn, coeff, p, int (idxOfMin), area);

//    Point * p_tangent = wgn->velocity->at (idxOfMin);

//    Point normal_out = {-p_tangent->y, p_tangent->x, 0};

//    *coeff = (normal_out | Point(*(wgn->position->at (idxOfMin)) - *p));

//    return int(idxOfMin);
//}

int Remainder (int dividend, int divisor)
{
    while (dividend >= divisor || dividend < 0)
        dividend += (dividend < 0 ? 1 : -1) * divisor;

    return dividend;
}

std::vector<std::vector <int>> Build (int i, std::pair<int, int> j, int k)
{
    int N = std::max(j.second, j.first) - std::min(j.second, j.first) +1;

    std::vector<std::vector <int>> out (3, std::vector<int>(size_t (N)));

    Eigen::ArrayXi x = Eigen::ArrayXi::LinSpaced (N, i, i);
    Eigen::ArrayXi y = Eigen::ArrayXi::LinSpaced (N, j.first, j.second);
    Eigen::ArrayXi z = Eigen::ArrayXi::LinSpaced (N, k, k);

    if (Remainder (k, 2) == 1)
        y.reverseInPlace ();

    out.at (0) = std::vector<int> (x.data(), x.data () + x.rows () * x.cols ());
    out.at (1) = std::vector<int> (y.data(), y.data () + y.rows () * y.cols ());
    out.at (2) = std::vector<int> (z.data(), z.data () + z.rows () * z.cols ());

    return out;
}

std::vector<std::vector <int>> Build (std::pair<int, int> i, int j, int k)
{
    int N = std::max(i.second, i.first) - std::min(i.second, i.first) +1;

    std::vector<std::vector <int>> out (3, std::vector<int>(size_t (N)));

    Eigen::ArrayXi x = Eigen::ArrayXi::LinSpaced (N, i.first, i.second);
    Eigen::ArrayXi y = Eigen::ArrayXi::LinSpaced (N, j, j);
    Eigen::ArrayXi z = Eigen::ArrayXi::LinSpaced (N, k, k);

    if (Remainder (k, 2) == 1)
        x.reverseInPlace ();

    out.at (0) = std::vector<int> (x.data(), x.data () + x.rows () * x.cols ());
    out.at (1) = std::vector<int> (y.data(), y.data () + y.rows () * y.cols ());
    out.at (2) = std::vector<int> (z.data(), z.data () + z.rows () * z.cols ());

    return out;
}
