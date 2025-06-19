#include <cstdlib>
#include <iostream>

int main() {
    std::cout << "Avvio esecuzione dello script R...\n";

    // Comando completo con path assoluto allo script
    const char* command =
        "/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe "
        "-Command \"& 'C:\\Program Files\\R\\R-4.5.0\\bin\\x64\\Rscript.exe' "
        "'C:\\Users\\francesca16\\TestTriangle\\Test_fmesher\\stella.R'\"";

    // Eseguo il comando, e reindirizzo output in log.txt per verifica
    std::string command_with_log = std::string(command) + " > r_output.txt 2>&1";

    std::cout << "Comando eseguito:\n" << command_with_log << "\n";

    int exit_code = std::system(command_with_log.c_str());

    std::cout << "Codice di uscita: " << exit_code << "\n";

    if (exit_code != 0) {
        std::cerr << "Errore durante l'esecuzione dello script R.\n";
        return 1;
    }

    std::cout << "Script R eseguito. Controlla r_output.txt per il log.\n";
    return 0;
}
