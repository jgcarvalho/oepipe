package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

var distributed bool

func conformers(ligandsDirPtr *string) {
	sdfs, _ := filepath.Glob(*ligandsDirPtr + "/*.sdf")
	fmt.Println(sdfs)
	os.Mkdir("./conformers", 0777)
	for _, f := range sdfs {
		outfile := "./conformers/" + filepath.Base(f)[:len(filepath.Base(f))-4] + ".oeb.gz"
		fmt.Println(outfile)
		// Change from 100 to 100000
		if distributed {
			cmd := exec.Command("srun", "omega2", "-in", f, "-out", outfile, "-sdEnergy", "-ewindow", "25.0", "-maxconfs", "100000", "-rms", "0.1")
			cmd.Start()
			defer cmd.Wait()
		} else {
			cmd := exec.Command("omega2", "-in", f, "-out", outfile, "-sdEnergy", "-ewindow", "25.0", "-maxconfs", "100", "-rms", "0.1")
			cmd.Run()
		}
	}
}

func docking(receptorsDirPtr *string) {
	conformers, _ := filepath.Glob("./conformers/*.oeb.gz")
	fmt.Println(conformers)
	outdir := "./docking/" + filepath.Base(*receptorsDirPtr) + "/"
	os.MkdirAll(outdir, 0777)
	for _, cf := range conformers {
		prefix := outdir + strings.TrimSuffix(filepath.Base(cf), ".oeb.gz")
		outfile := prefix + "_docked.oeb.gz"
		scorefile := prefix + "_score.txt"
		fmt.Println(outfile)
		if distributed {
			cmd := exec.Command("srun", "hybrid", "-receptor", *receptorsDirPtr+"/*", "-dbase", cf, "-docked_molecule_file", outfile, "-score_file", scorefile, "-dock_resolution", "High", "-num_poses", "25", "-save_component_scores", "-annotate_scores", "-prefix", prefix)
			cmd.Start()
			defer cmd.Wait()
		} else {
			// Change from Low to High
			cmd := exec.Command("hybrid", "-receptor", *receptorsDirPtr+"/*", "-dbase", cf, "-docked_molecule_file", outfile, "-score_file", scorefile, "-dock_resolution", "Low", "-num_poses", "25", "-save_component_scores", "-annotate_scores", "-prefix", prefix)
			cmd.Run()
		}
	}
}

func optimize(receptorsDirPtr *string) {
	docked, _ := filepath.Glob("./docking/" + filepath.Base(*receptorsDirPtr) + "/*_docked.oeb.gz")
	outdir := "./optimized/" + filepath.Base(*receptorsDirPtr) + "/"
	os.MkdirAll(outdir, 0777)
	var rec string
	for _, dck := range docked {
		rec = getBestReceptor(dck)
		prefix := outdir + strings.TrimSuffix(filepath.Base(dck), "_docked.oeb.gz")
		if distributed {
			cmd := exec.Command("srun", "szybki", "-p", rec, "-in", dck, "-prefix", prefix, "-residue", "1", "-protein_elec", "PB", "-am1bcc")
			cmd.Start()
			defer cmd.Wait()
		} else {
			// change residue from 1 to 3
			cmd := exec.Command("szybki", "-p", rec, "-in", dck, "-prefix", prefix, "-residue", "1", "-protein_elec", "PB", "-am1bcc")
			cmd.Run()
		}
	}
}

func entropyInSolution() {
	conformers, _ := filepath.Glob("./conformers/*.oeb.gz")
	outdir := "./entropy_insolution/"
	os.Mkdir(outdir, 0777)
	for _, cf := range conformers {
		prefix := outdir + strings.TrimSuffix(filepath.Base(cf), ".oeb.gz")
		fmt.Println(prefix)
		if distributed {
			cmd := exec.Command("srun", "szybki", "-entropy", "AN", "-sheffield", "-prefix", prefix, cf)
			cmd.Start()
			defer cmd.Wait()
		} else {
			cmd := exec.Command("szybki", "-entropy", "AN", "-sheffield", "-prefix", prefix, cf)
			cmd.Run()
		}
	}
}

func entropyBounded(receptorsDirPtr *string) {
	docked, _ := filepath.Glob("./docking/" + filepath.Base(*receptorsDirPtr) + "/*_docked.oeb.gz")
	outdir := "./entropy_bounded/" + filepath.Base(*receptorsDirPtr) + "/"
	os.MkdirAll(outdir, 0777)
	var rec string
	for _, dck := range docked {
		rec = getBestReceptor(dck)
		prefix := outdir + strings.TrimSuffix(filepath.Base(dck), "_docked.oeb.gz")
		if distributed {
			cmd := exec.Command("srun", "szybki", "-p", rec, "-entropy", "AN", "-prefix", prefix, dck)
			cmd.Start()
			defer cmd.Wait()
		} else {
			cmd := exec.Command("szybki", "-p", rec, "-entropy", "AN", "-prefix", prefix, dck)
			cmd.Run()
		}

	}

}

func getBestReceptor(dck string) string {
	scorefile := strings.Replace(dck, "_docked.oeb.gz", "_score.txt", 2)
	dat, err := ioutil.ReadFile(scorefile)
	if err != nil {
		panic(err)
	}
	fields := strings.Split(string(dat), "\t")

	return strings.TrimSuffix(fields[len(fields)-1], "\n")
}

func main() {
	flag.BoolVar(&distributed, "distributed", false, "distributed execution?")
	ligandsDirPtr := flag.String("l", "../Project/ligands", "ligands directory")
	receptorsDirPtr := flag.String("r", "../Project/receptor/K", "receptors directory")
	conformersPtr := flag.Bool("conformers", false, "generate conformers?")
	dockingPtr := flag.Bool("docking", false, "dock conformers?")
	optimizePtr := flag.Bool("optimize", false, "optimize conformers?")
	entropyISPtr := flag.Bool("eis", false, "entropy in solution?")
	entropyBPtr := flag.Bool("eb", false, "entropy bounded?")

	flag.Parse()

	fmt.Println("Distributed execution?", distributed)
	fmt.Println("Ligands directory:", *ligandsDirPtr)
	fmt.Println("Receptors directory:", *receptorsDirPtr)
	fmt.Println("Generate conformers?", *conformersPtr)
	fmt.Println("Dock conformers?", *dockingPtr)
	fmt.Println("Optimize docked conformers", *optimizePtr)

	if *conformersPtr {
		fmt.Println("Generating conformers")
		conformers(ligandsDirPtr)
	}
	if *dockingPtr {
		fmt.Println("Docking ligands")
		docking(receptorsDirPtr)
	}
	if *entropyISPtr {
		fmt.Println("Estimating entropy in solution")
		entropyInSolution()
	}
	if *entropyBPtr {
		fmt.Println("Estimating entropy bounded")
		entropyBounded(receptorsDirPtr)
	}
	if *optimizePtr {
		optimize(receptorsDirPtr)
	}

}