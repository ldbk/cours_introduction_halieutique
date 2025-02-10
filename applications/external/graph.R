# Installer et charger ggplot2 si ce n'est pas déjà fait
# install.packages("ggplot2")
library(ggplot2)

# Création de données fictives
annees <- 2000:2020
# Exemple d'évolution de la biomasse (B)
B <- seq(1.0, 1.5, length.out = length(annees))
# Blim et Bmsy définis comme des constantes (ou évoluant très peu)
Blim <- rep(1.2, length(annees))
Bmsy <- rep(1.3, length(annees))
# Exemple d'évolution de la mortalité par pêche (F)
F <- seq(0.2, 0.1, length.out = length(annees))

# Rassembler les données dans un data.frame
data <- data.frame(annee = annees, B = B, Blim = Blim, Bmsy = Bmsy, F = F)

# Pour représenter F sur le même graphique, on définit un facteur de conversion
# afin de mettre F à l'échelle des valeurs de B.
facteur <- max(data$B) / max(data$F)

# Création du graphique avec ggplot2
p <- ggplot(data, aes(x = annee)) +
  # Courbe de la biomasse B
  geom_line(aes(y = B, color = "B"), size = 1) +
  # Courbe de Blim (en tirets)
  geom_line(aes(y = Blim, color = "Blim"), linetype = "dashed", size = 1) +
  # Courbe de Bmsy (en pointillés)
  geom_line(aes(y = Bmsy, color = "Bmsy"), linetype = "dotted", size = 1) +
  # Courbe de la mortalité par pêche F (mise à l'échelle)
  geom_line(aes(y = F * facteur, color = "Mortalité par pêche"), linetype = "dotdash", size = 1) +
  # Définition des axes : l'axe principal pour la biomasse et un axe secondaire pour F
  scale_y_continuous(name = "Biomasse",
                     sec.axis = sec_axis(~ . / facteur, name = "Mortalité par pêche")) +
  labs(x = "Année", color = "Paramètres") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Affichage du graphique
print(p)
