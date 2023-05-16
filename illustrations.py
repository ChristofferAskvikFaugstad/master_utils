from utils.master import *


def add_path(
    diagram: ED,
    names: List[str],
    relative_energies: List[float],
    energy_barriers: List[float],
    state_seperations: List[float],
    start_height: float,
    start_position: float,
    start_link: int,
    color,
    include_start: bool = False,
    top_relative = True
):
    height = start_height
    position = start_position
    TScolor = color
    org_num_links = len(diagram.energies)

    if include_start:
        pass
    else:
        for i, (name, dE, Eb) in enumerate(
            zip(names[1:], relative_energies, energy_barriers)
        ):
            position += state_seperations[2 * i]
            diagram.add_level(
                Eb + height,
                bottom_text="",
                linestyle="-",
                color=TScolor,
                position=position,
                top_text= f"{Eb + height}" if top_relative else rf"{Eb}$^\ddag$",
                )
            position += state_seperations[2 * i + 1]
            height += dE
            diagram.add_level(
                height,
                bottom_text=name,
                linestyle="-",
                color=TScolor,
                position=position,
            )
        for i in range(org_num_links, len(diagram.energies) - 1):
            diagram.add_link(i, i + 1, color=color)
        diagram.add_link(start_link, org_num_links, color=color)


def add_CO_activation(
    ax,
    CO_Ead,
    colors,
    namess,
    relative_energiess,
    energy_barrierss,
    state_seperationss,
    title,
    labels=["CO", "COH", "HCO", "HCOH"],
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    loc="upper center",
    diagram_ratio=3,
    diagram_offset_ratio=0.01,
    diagram_round_energies_at_digit=3,
    top_relative=True,
):
    diagram = ED()
    diagram.offset_ratio = diagram_offset_ratio
    diagram.ratio = diagram_ratio
    diagram.round_energies_at_digit = diagram_round_energies_at_digit
    handels = []
    diagram.add_level(0.0, "CO", color="black")
    diagram.add_level(CO_Ead, "CO*", color="black")
    diagram.add_link(0, 1, color="black")
    for (
        names,
        relative_energies,
        energy_barriers,
        state_seperations,
        color,
        label,
    ) in zip(
        namess, relative_energiess, energy_barrierss, state_seperationss, colors, labels
    ):
        add_path(
            diagram,
            names=names,
            relative_energies=relative_energies,
            energy_barriers=energy_barriers,
            state_seperations=state_seperations,
            start_height=CO_Ead,
            start_position=2.0,
            start_link=1,
            color=color,
            top_relative=top_relative,
        )
        handels.append(plt.Line2D([], [], color=color, marker=" ", ls="-", label=label))
    diagram.plot(ylabel="Relative Energy [eV]", ax=ax)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax.legend(
        handles=handels,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        fancybox=True,
        shadow=True,
        ncol=ncol,
    )
    ax.set_title(title, pad=20)


def add_CO_activation_by_path(
    ax,
    CO_Ead,
    colors,
    namess,
    relative_energiess,
    energy_barrierss,
    state_seperationss,
    title,
    labels=["CO", "COH", "HCO", "HCOH"],
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    loc="upper center",
    diagram_ratio=3,
    diagram_offset_ratio=0.01,
    diagram_round_energies_at_digit=3,
    top_relative=True,
):
    diagram = ED()
    diagram.offset_ratio = diagram_offset_ratio
    diagram.ratio = diagram_ratio
    diagram.round_energies_at_digit = diagram_round_energies_at_digit
    handels = []
    diagram.add_level(0.0, "CO", position=0, color="black")
    for (
        names,
        relative_energies,
        energy_barriers,
        state_seperations,
        color,
        label,
    ) in zip(
        namess, relative_energiess, energy_barrierss, state_seperationss, colors, labels
    ):
        i = len(diagram.positions)
        diagram.add_level(
            relative_energies[0], names[0], position=state_seperations[0], color=color
        )
        diagram.add_link(0, i, color=color)
        add_path(
            diagram,
            names=names,
            relative_energies=relative_energies[1:],
            energy_barriers=energy_barriers,
            state_seperations=state_seperations[1:],
            start_height=CO_Ead,
            start_position=state_seperations[0],
            start_link=i,
            color=color,
        )
    diagram.plot(ylabel="Relative Energy [eV]", ax=ax)
    ax.set_title(title, pad=20)
