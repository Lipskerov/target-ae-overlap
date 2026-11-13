# Prints the openFDA API key from the FDA_API_KEY environment variable (empty if unset; openFDA then applies the
# keyless daily limit). No key is stored in this repository.
import os
print(os.environ.get("FDA_API_KEY", ""))
