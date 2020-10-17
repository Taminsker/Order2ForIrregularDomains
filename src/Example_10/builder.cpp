#include "headers.h"

// Utile

double phi (Point p, double t)
{
//    (void)t;

    //    return fabs(p.x - 0.5) - 0.1;
//    return fabs(p.x - 0.5 - std::exp (-p.x + 0.5) * (std::exp(t) - 1)) - 0.1;

//    double i = std::cos(5. * M_PI * (p.x - t)) - std::cos(5. * M_PI * p.x);
    p.x -= u(p, 0) - u(p, t);
    return fabs(p.x - 0.5) - 0.1;

}

double f (Point a, double t)
{
    (void) t;
    (void) a;
    return 0;
}

double u (Point a, double t)
{
//    return std::exp(t - a.x + 0.5) - 1.;
    return std::sin(5. * M_PI * (a.x - t));
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

void Make(Point* minima, Point* extrema, double time_final, STYLE style, std::vector<int> listNx, std::vector<int> listNy, std::vector<int> listNz, double h0)
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

//        std::cout.clear();
        mesh->Build ();
        mesh->Print ();
//        std::cout.setstate(std::ios_base::failbit);


        double hx = mesh->Get_hx ();
        double hy = mesh->Get_hy ();
        double hz = mesh->Get_hz ();

        if (hx < 1e-20) {hx = 1.;}
        if (hy < 1e-20) {hy = 1.;}
        if (hz < 1e-20) {hz = 1.;}

        switch (style)
        {
        case BackwardEuler1:
            dt = 1. * std::min(hx, std::min (hy, hz));
            name = "BackwardEuler1";
            break;
        case BackwardEuler2:
            dt = std::min(hx * hx, std::min (hy * hy, hz * hz));
            name = "BackwardEuler2";
            break;
        case CrankNicolson:
            dt = 1. * std::min(hx, std::min (hy, hz));
            name = "CrankNicolson";
            break;
        }

        int IterMax = int (time_final / dt);

        listdt.push_back (dt);

        // Construction des vecteurs phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi, 0.);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);
        phi_vec = FunToVec (mesh, phi, 0.);

        mesh->Print ();

        // u_num
        Vector u_num = FunToVec (mesh, u, 0.);

        // Field;
        Field* field = GetWField (mesh, &phi_vec, &u_num, &listPoint, h0);
        Vector NormGradPhi = ComputeNorm(&field->GradPhi);

        // Écriture dans des fichiers
        Writer writer (mesh);
        writer.SetFilename (std::string("Example_10_") + name + std::string ("_Nx_") + std::to_string (Nx) + std::string ("_Ny_") + std::to_string (Ny) + std::string ("_Nz_") + std::to_string (Nz));

        bool mybool = true;

        for (int it = 0; (it < IterMax) && mybool; ++it)
        {
            double time = double(it) * dt;

            std::cout.clear();
            std::cout << "\r" << INDENT << "Iteration : time : " << time << "/" << time_final << " dt : " << dt << "            " << std::flush;
            std::cout.setstate(std::ios_base::failbit);

            if (time > 0)
            {
                Matrix Al, Ar;
                GetMatrix (mesh, dt, style, &Al, &Ar);

                Extrapole (mesh, &u_num);

                // Matrix
                Matrix A = Al;
                Vector b = Ar * u_num;

                ImposeDirichlet (mesh, &A, &b, u, listPoint, time);

                u_num = Solve (A, b, IMPLICIT);

            }

            Vector u_ana = FunToVec (mesh, u, time);

            mesh->MakeZeroOnExternOmegaInVector (&u_ana);
            mesh->MakeZeroOnExternOmegaInVector (&u_num);

            Vector err_abs = GetErrorAbs (mesh, u_ana, u_num);

            if (IterMax-1 == it)
            {
                err_l1.push_back (GetErrorl1 (mesh, u_ana, u_num)); // erreur l1
                err_linf.push_back (GetErrorlinf (mesh, u_ana, u_num)); // erreur linf
                err_rela.push_back (GetErrorRela (mesh, u_ana, u_num));// erreur relative
                Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
                h.push_back (std::sqrt(p|p));
            }

            // Extrapole partie
            Extrapole (mesh, &NormGradPhi);
            Extrapole (mesh, &field->W);
            Extrapole (mesh, &field->Normals);
            Extrapole (mesh, &field->GradPhi);
            Extrapole (mesh, &field->GradTemperature);

            // Écriture
            writer.SetCurrentIteration (it);
            writer.SetVectorNumerical (&u_num);
            writer.SetVectorAnalytical (&u_ana);
            writer.SetVectorErrorAbs (&err_abs);
            writer.SetVectorPhi (&phi_vec);
            writer.SetNormPhi (&NormGradPhi);
            writer.SetVectorW_new (&field->W);
            writer.SetVectorNormals (&field->Normals);
            writer.SetVectorGradPhi (&field->GradPhi);
            writer.SetVectorGradTemperature (&field->GradTemperature);
            writer.WriteNow ();


            //            // Compute for Next iteration
            //            delete field;
            //            listPoint = mesh->GetListOfIndexPoints ();

            //            field = GetWField (mesh, &phi_vec, &u_num, &listPoint, h0);

            //            //        dt = Compute_dt (mesh, field);

            //            NormGradPhi = ComputeNorm(&field->GradPhi);

            //            phi_vec = IteratePhi (mesh, &field->W, dt, &phi_vec);

            phi_vec = FunToVec (mesh, phi, time);
            mesh->RemoveAllNotCartesianPoints ();
            listPoint = MakeBorderPoints (mesh, &phi_vec);
            phi_vec = FunToVec (mesh, phi, time);

            //            std::cout.setstate(std::ios_base::failbit);

            //            phi_vec = ReInitPhi (mesh, &phi_vec);

            //            Extrapole (mesh, &phi_vec);
        }


        delete field;
        delete mesh;
    }

    std::cout.clear();

    std::cout << "\r#Summary " << name << " (time final = " << time_final  << ")            " << std::endl;

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
