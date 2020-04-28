#include "headers.h"

// Utile

double phi (Point p, double t)
{
    (void)t;

    return EuclidianDist (p, Point(0.5, 0.5, 0.5)) - 0.15;
}

double f (Point a, double t)
{
    (void) t;
    (void) a;
    return 0;
}

double u (Point a, double t)
{
    return std::exp(- 3. * t) * std::sin(a.x) * std::sin(a.y) * std::sin(a.z);
}

// Fonctions builder

void GetMatrix(Mesh* mesh, double dt, STYLE style, Matrix* Al, Matrix* Ar)
{
    int N = mesh->GetNumberOfTotalPoints ();

    Matrix A = - Laplacian (mesh);
    RemovePeriodicity (mesh, &A);
    Matrix Id (N, N);
    Id.setIdentity ();

    switch (style)
    {
    case BackwardEuler1:
    case BackwardEuler2:
        *Al = (1./ dt) * Id + A;
        *Ar = (1./ dt) * Id;
        return;
    case CrankNicolson:
        *Al = (1./ dt) * Id + 0.5 * A;
        *Ar = (1./ dt) * Id - 0.5 * A;
        return;
    }
}

void Make(Point* minima, Point* extrema, int Nt, STYLE style, std::vector<int> listNx, std::vector<int> listNy, std::vector<int> listNz)
{

    std::string name = "";

    double dt = 0;

    if (    listNx.size () != listNy.size () ||
            listNx.size () != listNz.size () ||
            listNy.size () != listNz.size ())
        return;

    std::vector<double> err_l1 = {}; // erreur l1
    std::vector<double> err_linf = {}; // erreur linf
    std::vector<double> err_rela = {}; // erreur relative
    std::vector<double> h = {}; // pas d'espace h
    std::vector<double> listdt = {}; // pas d'espace h


    size_t n = listNx.size();

    std::cout.setstate(std::ios_base::failbit);

    for (size_t i = 0; i < n; ++i)
    {
        int Nx = listNx.at (i);
        int Ny = listNy.at (i);
        int Nz = listNz.at (i);

        Mesh* mesh = new Mesh();
        mesh->SetBounds (new Point(*minima), new Point(*extrema));
        mesh->Set_Nx (std::max(Nx, 1));
        mesh->Set_Ny (std::max(Ny, 1));
        mesh->Set_Nz (std::max(Nz, 1));

        mesh->Build ();
        mesh->Print ();

        double hx = mesh->Get_hx ();
        double hy = mesh->Get_hy ();
        double hz = mesh->Get_hz ();

        if (hx < 1e-20) {hx = 1.;}
        if (hy < 1e-20) {hy = 1.;}
        if (hz < 1e-20) {hz = 1.;}

        switch (style)
        {
        case BackwardEuler1:
            dt = std::min(hx, std::min (hy, hz));
            name = "BackwardEuler1";
            break;
        case BackwardEuler2:
            dt = std::min(hx * hx, std::min (hy * hy, hz * hz));
            name = "BackwardEuler2";
            break;
        case CrankNicolson:
            dt = std::min(hx, std::min (hy, hz));
            name = "CrankNicolson";
            break;
        }

        listdt.push_back (dt);

        // Construction de vecteur phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);

        Matrix Al;
        Matrix Ar;

        GetMatrix (mesh, dt, style, &Al, &Ar);

        Vector u_num, u_ana;
        u_num = u_ana = FunToVec (mesh, u, 0.);

        Vector err_abs = GetErrorAbs (mesh, u_ana, u_num);

        int compt = 0;

        Writer writer (mesh);
        writer.SetFilename (std::string("Example_9_") + name + std::string ("_Nx_") + std::to_string (Nx) + std::string ("_Ny_") + std::to_string (Ny) + std::string ("_Nz_") + std::to_string (Nz));

        writer.SetCurrentIteration (compt); // Itérations lorsqu'il y a du temps
        writer.SetVectorNumerical (&u_num);
        writer.SetVectorAnalytical (&u_ana);
        writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
        writer.SetVectorErrorAbs (&err_abs);

        writer.WriteNow ();
        compt++;


        for (int it = 1; it <= Nt; it++)
        {
            std::cout.clear();
            std::cout << "\r" << INDENT << "Iteration : index : " << i << "/" << n << " it : " << it << "/" << Nt << "                 "<< std::flush;
            std::cout.setstate(std::ios_base::failbit);

            double t = double(it) * dt;

            Matrix A = Al;
            Vector b = Ar * u_num;

            ImposeDirichlet (mesh, &A, &b, u, listPoint, t);

            u_num = Solve (A, b, IMPLICIT);
            u_ana = FunToVec (mesh, u, t);

            mesh->MakeZeroOnExternOmegaInVector (&u_ana);
            mesh->MakeZeroOnExternOmegaInVector (&u_num);

            err_abs = GetErrorAbs (mesh, u_ana, u_num);

            if (it == Nt - 1)
            {
                err_l1.push_back (GetErrorl1 (mesh, u_ana, u_num));
                err_linf.push_back (GetErrorlinf (mesh, u_ana, u_num));
                err_rela.push_back (GetErrorRela (mesh, u_ana, u_num));

                Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
                h.push_back (std::sqrt(p|p));
            }



            if (it%1 == 0)
            {
                //                 Écriture dans des fichiers
                writer.SetCurrentIteration (compt); // Itérations lorsqu'il y a du temps
                writer.SetVectorNumerical (&u_num);
                writer.SetVectorAnalytical (&u_ana);
                writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
                writer.SetVectorErrorAbs (&err_abs);

                writer.WriteNow ();
                compt++;
            }
        }


        delete mesh;
    }

    std::cout.clear();

    std::cout << "\r#Summary " << name << " (Nt = " << Nt-1 << ")            " << std::endl;

    std::cout << "Nx            : " << listNx << std::endl;
    std::cout << "Ny            : " << listNy << std::endl;
    std::cout << "Nz            : " << listNz << std::endl;
    std::cout << "dt            : " << listdt << std::endl;
    std::cout << "l1-error      : " << err_l1 << std::endl;
    std::cout << "Order         : " << Order(err_l1, h) << std::endl;
    std::cout << "linf-error    : " << err_linf << std::endl;
    std::cout << "Order         : " << Order(err_linf, h) << std::endl;
    std::cout << "rela-error    : " << err_rela << std::endl;
    std::cout << "Order         : " << Order(err_rela, h) << std::endl;

    std::cout << std::endl;
    return;
}
